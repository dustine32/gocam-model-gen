import gocamgen
from gocamgen.gocamgen import expand_uri_wrapper, ACTS_UPSTREAM_OF_RELATIONS
import unittest
import logging
from gocamgen.filter_rule import WBFilterRule, MGIFilterRule
from gen_models_by_gene import AssocExtractor, GoCamBuilder
from gocamgen.triple_pattern_finder import TriplePattern, TriplePatternFinder, TriplePair, TriplePairCollection
from rdflib.term import URIRef
from gocamgen.rdflib_sparql_wrapper import RdflibSparqlWrapper
from gocamgen.subgraphs import AnnotationSubgraph
from gocamgen.utils import ShexHelper

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("gocamgen.gocamgen")
logger.setLevel(logging.ERROR)


class TestReferencePreference(unittest.TestCase):

    def test_ref_picker(self):
        test_refs = [
            "GO_REF:0000483",
            "doi:485930",
            "WB_REF:WBPaper00003384",
            "PMID:9834189",
        ]
        ref_picker = gocamgen.gocamgen.ReferencePreference()

        result = ref_picker.pick(test_refs)
        self.assertEqual(result, "PMID:9834189")


class TestGoCamModel(unittest.TestCase):
    BUILDER = GoCamBuilder()  # Takes a sec to init so only make once

    def setUp(self):
        pass

    def gen_model(self, gpad_file, test_gene, filter_rule):
        extractor = AssocExtractor(gpad_file, filter_rule)
        assocs_by_gene = extractor.group_assocs()

        if test_gene not in assocs_by_gene:
            self.fail("ERROR: specific gene {} not found in filtered annotation list".format(test_gene))
        else:
            model = TestGoCamModel.BUILDER.translate_to_model(test_gene, assocs_by_gene[test_gene])
            if model:
                return model
            else:
                self.fail("Couldn't generate model for {}".format(test_gene))

    def test_triple_finder(self):
        test_gene = "WB:WBGene00006498"
        model = self.gen_model(gpad_file="resources/test/wb_6498.gpad", test_gene=test_gene,
                          filter_rule=WBFilterRule())
        # Get model.writer.graph whatever and check for loose evidence (not attached to axioms)
        # Orphaned evidence - how to I find these in debugger? They're in rdflib writer somewhere

        # Can this work?
        #   (MF, ENABLED_BY, GP) & (Same MF, has input, GP)
        pattern_a = TriplePattern([("GO:0003674", gocamgen.gocamgen.ENABLED_BY, test_gene)])
        pattern_b = TriplePattern([("GO:0003674", URIRef(expand_uri_wrapper("BFO:0000050")), "GO:0019953")])
        whole_pattern = TriplePattern([("GO:0003674", gocamgen.gocamgen.ENABLED_BY, test_gene),
                                 ("GO:0003674", URIRef(expand_uri_wrapper("BFO:0000050")), "GO:0019953")])

        triple_finder = TriplePatternFinder()
        a_triples = triple_finder.find_pattern_recursive(model, pattern_a)
        # print("A count: {}".format(len(a_triples)))
        b_triples = triple_finder.find_pattern_recursive(model, pattern_b)
        # print("B count: {}".format(len(b_triples)))
        found_chains = triple_finder.find_pattern_recursive(model, whole_pattern)
        # print(found_chains)
        # print("Chain count: {}".format(len(found_chains)))
        # for fc in found_chains:
        #     print(contract_uri_wrapper(model.individual_label_for_uri(fc[1][2])[0])[0])

        triple_pair = TriplePair(pattern_a.ordered_triples[0], pattern_b.ordered_triples[0],
                                 connecting_entity="GO:0003674")
        tp_collection = TriplePairCollection()
        tp_collection.chain_collection.append(triple_pair)
        uri_tp_collection = triple_finder.find_connected_pattern(model, tp_collection)

        self.assertGreaterEqual(len(uri_tp_collection.chain_collection), 1)

    def test_evidence(self):
        # gpad_file = "wb_903.gpad"
        # test_gene = "WB:WBGene00000903"
        # filter_rule = WBFilterRule()

        self.assertEqual(1, 1)

    def test_has_input(self):
        # See https://github.com/geneontology/gocamgen/issues/39#issuecomment-479988904 for background
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())
        # Look for translation of 'GO:0000977 has_direct_input(WB:WBGene00036254)'
        found_triples = model.triples_by_ids("GO:0000977", URIRef(expand_uri_wrapper("RO:0002233")),
                                             "WB:WBGene00036254")
        self.assertGreaterEqual(len(found_triples), 1, "No has_input extensions translated")

    def test_extension_pipe_separation(self):
        # See https://github.com/geneontology/gocamgen/issues/40
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())

        # Look for count of 'WB:WBGene00003167 contributes_to GO:0000977'
        found_triples = model.triples_by_ids("WB:WBGene00003167", gocamgen.gocamgen.CONTRIBUTES_TO,
                                             "GO:0000977")
        self.assertGreaterEqual(len(found_triples), 3,
                                "Less than 3 annotations for WB:WBGene00003167 contributes_to GO:0000977")

    def test_no_dup_individuals(self):
        # See https://github.com/geneontology/gocamgen/issues/40
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_2159711", test_gene="MGI:MGI:2159711",
                          filter_rule=MGIFilterRule())

        # Look for 'MGI:MGI:2159711 LOCATED_IN GO:0044297'. Should only be one.
        found_triples = model.triples_by_ids("MGI:MGI:2159711", gocamgen.gocamgen.LOCATED_IN,
                                             "GO:0044297")
        self.assertEqual(len(found_triples), 1)

        # Case of multiple experimental annotations to same GP and term.
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_98956", test_gene="MGI:MGI:98956",
                               filter_rule=MGIFilterRule())

        # This example has like 10-20 exp assertions for MGI:MGI:98956 to GO:0060070
        assertion_count = len(model.associations.collapsed_associations)
        sparql_wrapper = RdflibSparqlWrapper()
        gp = "MGI:MGI:98956"
        term = "GO:0060070"
        qres = sparql_wrapper.find_involved_in_translated(model.graph, gp, term)
        self.assertEqual(len(qres), assertion_count)

    def test_has_regulation_target(self):
        # Examples:
        # F - MGI:MGI:107771 GO:0005096 'has_regulation_target(MGI:MGI:97846)|has_regulation_target(MGI:MGI:2180784)'
        # P - WB:WBGene00013591 GO:0042594 'causally_upstream_of(GO:0001934),has_regulation_target(WB:WBGene00008480)'
        # Which has_regulation_target bucket does this fall into? None so far (GO:0042594 is "response to starvation")
        # bucket = gocamgen.gocamgen.has_regulation_target_bucket(ont, "GO:0001934")

        # model = gen_model(gpad_file="resources/test/wb.gpad.WBGene00013591", test_gene="WB:WBGene00013591",
        #                   filter_rule=WBFilterRule())
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())

        self.assertEqual(1, 1)

    def test_acts_upstream_of(self):
        # TODO: Test for MGI:MGI:1206591
        self.assertEqual(1, 1)

        # what about MGI:MGI:1914305?
        # resources/test/mgi.gpa.MGI_1914305
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_1914305", test_gene="MGI:MGI:1914305",
                               filter_rule=MGIFilterRule())

        sparql_wrapper = RdflibSparqlWrapper()
        gp = "MGI:MGI:1914305"
        term = "GO:0007416"
        causally_relation = model.get_causally_upstream_relation(ACTS_UPSTREAM_OF_RELATIONS["acts_upstream_of_or_within"])
        qres = sparql_wrapper.find_acts_upstream_of_translated(model.graph, gp, causally_relation, term)
        self.assertEqual(len(qres), 1)

    def test_sparql(self):
        # Just gonna see what we can sparql out of this guy.
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                               filter_rule=WBFilterRule())

        sparql_wrapper = RdflibSparqlWrapper()
        gp = "WB:WBGene00003167"
        term = "GO:0007638"
        qres = sparql_wrapper.find_involved_in_translated(model.graph, gp, term)

        self.assertEqual(len(qres), 1)

    def test_graph_to_sparql(self):
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                               filter_rule=WBFilterRule())

        # Translation of a regulation of gene expression has_regulation_target extension
        # Write this:
        # GO:0003674-1 enabled_by (RO:0002333) WB:WBGene00003167-1
        # GO:0003674-1 part_of (BFO:0000050) GO:0045944-1
        # GO:0003674-1 has_input (RO:0002233) WB:WBGene00003167-2
        # GO:0003674-2 enabled_by (RO:0002333) WB:WBGene00003167-2
        # GO:0003674-1 causally upstream of, positive effect (RO:0002304) GO:0003674-2
        g = AnnotationSubgraph({"source_line": "fake annotation"})
        mf_a = g.add_instance_of_class("GO:0003674")
        gp_a = g.add_instance_of_class("WB:WBGene00003167")
        bp_a = g.add_instance_of_class("GO:0045944")
        gp_b = g.add_instance_of_class("WB:WBGene00003167")
        mf_b = g.add_instance_of_class("GO:0003674")

        g.add_edge(mf_a, "RO:0002333", gp_a)
        g.add_edge(mf_a, "BFO:0000050", bp_a)
        g.add_edge(mf_a, "RO:0002233", gp_b)
        g.add_edge(mf_b, "RO:0002333", gp_b)
        g.add_edge(mf_a, "RO:0002304", mf_b)

        # g.print_matches_in_model(model)
        results = g.find_matches_in_model(model)

        # self.assertEqual(len(results), 1)
        self.assertEqual(1, 1)  # TODO: This example extension is going to be converted to has_input in near future

    def test_occurs_in_nesting(self):
        # Generate example of GO:CC->CL->EMAPA nesting for MF or BP primary terms
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_1336882_occurs_in", test_gene="MGI:MGI:1336882",
                               filter_rule=MGIFilterRule())
        # Look for chain of BP-occurs_in->CC-part_of->CL-part_of->EMAPA
        sparql_wrapper = RdflibSparqlWrapper()
        res = sparql_wrapper.find_nested_location_chain(model.graph, "BiologicalProcess", "GO:0070625", "GO:0045178", "CL:0002064", "EMAPA:35651")
        self.assertEqual(len(res), 1, "No nested chain found in MGI:MGI:1336882 test model")

        # Split example annotation with extension set of occurs_in(CL),occurs_in(CL),occurs_in(EMAPA) and nest
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_1915585", test_gene="MGI:MGI:1915585",
                               filter_rule=MGIFilterRule())
        # Look for 2chainz: BP-occurs_in->CL:0000589-part_of->EMAPA:17597 and
        #                   BP-occurs_in->CL:0000601-part_of->EMAPA:17597
        res = sparql_wrapper.find_nested_location_chain(model.graph, "BiologicalProcess", "GO:0090102", "CL:0000589", "EMAPA:17597")
        self.assertEqual(len(res), 1)
        res = sparql_wrapper.find_nested_location_chain(model.graph, "BiologicalProcess", "GO:0090102", "CL:0000601", "EMAPA:17597")
        self.assertEqual(len(res), 1)

    def test_part_of_nesting(self):
        # Generate example of GO:CC->CL->EMAPA nesting for CC primary terms
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_1336882_part_ofs", test_gene="MGI:MGI:1336882",
                               filter_rule=MGIFilterRule())
        sparql_wrapper = RdflibSparqlWrapper()
        # Look for chain of primaryCC-part_of->CC-part_of->CL-part_of->EMAPA
        res = sparql_wrapper.find_nested_location_chain(model.graph, "CellularComponent", "GO:0042588", "GO:1990794", "CL:0002064",
                                                        "EMAPA:35651")
        self.assertEqual(len(res), 1, "No nested CC-part_of-> chain found in MGI:MGI:1336882 test model")

        # This example has CC primary term SNARE complex, which should affect the GP->primaryCC->CC relations
        res = sparql_wrapper.find_nested_location_chain(model.graph, "ProteinContainingComplex", "GO:0031201", "GO:0042589", "CL:0002064", "EMAPA:35651")
        self.assertEqual(len(res), 1, "No nested complex-located_in->part_of-> chain found in MGI:MGI:1336882 test model")

    def test_with_from_evidence(self):
        # Check to ensure with/from values are appropriately added to evidence or translated to has_input edges
        sparql_wrapper = RdflibSparqlWrapper()
        enabled_by = "RO:0002333"
        has_input = "RO:0002233"

        # Non-protein binding case - should have evidence-with property on evidence individual:
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.mgi_1859682", test_gene="MGI:MGI:1859682",
                               filter_rule=MGIFilterRule())
        qres = sparql_wrapper.find_evidence_with(model.graph, "GO:0008520", enabled_by, "MGI:MGI:1859682")
        withs = []
        for r in qres:
            withs.append(str(r["evi_with"]))
        self.assertIn("MGI:MGI:2180333", withs)

        # Non-protein binding case w/ commas
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_87859", test_gene="MGI:MGI:87859",
                               filter_rule=MGIFilterRule())
        qres = sparql_wrapper.find_evidence_with(model.graph, "GO:0004713", enabled_by, "MGI:MGI:87859")
        withs = []
        for r in qres:
            withs.append(str(r["evi_with"]))
        # Current test is invalid now that we aren't splitting commas.
        #  AssertionError: 'MGI:MGI:3525886' not found in ['MGI:MGI:87860', 'MGI:MGI:87860', 'MGI:MGI:2176222,MGI:MGI:3525886', 'MGI:MGI:104738', 'MGI:MGI:87860']
        # self.assertIn("MGI:MGI:3525886", withs)
        # self.assertIn("MGI:MGI:2176222", withs)

        # Protein binding case - should have with/from value as term-has_input->with/from:
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00000018", test_gene="WB:WBGene00000018",
                               filter_rule=WBFilterRule())
        qres = sparql_wrapper.find_triple_by_class(model.graph, "GO:0005515", has_input, "WB:WBGene00015146")
        self.assertEqual(len(qres), 1)

        # DNA binding case - same as protein binding - should have with/from value as term-has_input->with/from:
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_97747", test_gene="MGI:MGI:97747",
                               filter_rule=MGIFilterRule())
        qres = sparql_wrapper.find_triple_by_class(model.graph, "GO:0003677", has_input, "UniProtKB:P28700")
        self.assertEqual(len(qres), 1)

        # Protein binding descendant case - should have with/from value as term-has_input->with/from:
        #  GO:0019894 - kinesin binding
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00000099", test_gene="WB:WBGene00000099",
                               filter_rule=WBFilterRule())
        qres = sparql_wrapper.find_triple_by_class(model.graph, "GO:0019894", has_input, "WB:WBGene00006974")
        self.assertEqual(len(qres), 1)

        # Protein binding w/ pipes
        #  GO:0017022 - myosin binding
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00002173", test_gene="WB:WBGene00002173",
                               filter_rule=WBFilterRule())
        qres = sparql_wrapper.find_triple_by_class(model.graph, "GO:0017022", has_input)
        # Assert there are 3 distinct individuals of GO:0017022 in these results
        result_subj_individuals = []
        for r in qres:
            subject_iri = r["s"]
            if subject_iri not in result_subj_individuals:
                result_subj_individuals.append(subject_iri)
        self.assertEqual(len(result_subj_individuals), 3)

    def test_iso_mgi_filtering(self):
        extractor = AssocExtractor(gpad_file="resources/test/mgi.gpa.MGI_1920081", filter_rule=MGIFilterRule())
        filtered_assocs = []
        for a in extractor.assocs:
            if extractor.assoc_filter.validate_line(a):
                filtered_assocs.append(a)
        self.assertEqual(len(filtered_assocs), 1)

    def test_complex_shex_shape_recognition(self):
        complex_root = "GO:0032991"
        shex_helper = ShexHelper()

        shape = shex_helper.shape_from_class(complex_root, TestGoCamModel.BUILDER.ext_mapper.go_aspector)
        self.assertEqual(shape, "ProteinContainingComplex")

        shape = shex_helper.shape_from_class("GO:0003674", TestGoCamModel.BUILDER.ext_mapper.go_aspector)
        self.assertEqual(shape, "MolecularFunction")


if __name__ == '__main__':
    unittest.main()