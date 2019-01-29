from ontobio.rdfgen.assoc_rdfgen import CamRdfTransform, TurtleRdfWriter, genid, prefix_context
from ontobio.vocabulary.relations import OboRO, Evidence
from ontobio.vocabulary.upper import UpperLevel
# from ontobio.util.go_utils import GoAspector
from prefixcommons.curie_util import expand_uri
from rdflib.namespace import OWL, RDF
from rdflib import Literal
from rdflib.term import URIRef
from rdflib.namespace import Namespace
import rdflib
import networkx
import logging
import argparse
import datetime
import os.path as path

# logging.basicConfig(level=logging.INFO)

ro = OboRO()
evt = Evidence()
upt = UpperLevel()
LEGO = Namespace("http://geneontology.org/lego/")
LAYOUT = Namespace("http://geneontology.org/lego/hint/layout/")
PAV = Namespace('http://purl.org/pav/')
DC = Namespace("http://purl.org/dc/elements/1.1/")

### ontobio/dustine32-issue-203 in here - delete once https://github.com/biolink/ontobio/pull/281 is in ontobio ###
def get_ancestors(ontology, go_term):
    all_ancestors = ontology.ancestors(go_term)
    all_ancestors.append(go_term)
    subont = ontology.subontology(all_ancestors)
    return subont.ancestors(go_term, relations=["subClassOf","BFO:0000050"])

def is_biological_process(ontology, go_term):
    bp_root = "GO:0008150"
    if go_term == bp_root:
        return True
    ancestors = get_ancestors(ontology, go_term)
    if bp_root in ancestors:
        return True
    else:
        return False

def is_molecular_function(ontology, go_term):
    mf_root = "GO:0003674"
    if go_term == mf_root:
        return True
    ancestors = get_ancestors(ontology, go_term)
    if mf_root in ancestors:
        return True
    else:
        return False

def is_cellular_component(ontology, go_term):
    cc_root = "GO:0005575"
    if go_term == cc_root:
        return True
    ancestors = get_ancestors(ontology, go_term)
    if cc_root in ancestors:
        return True
    else:
        return False

def go_aspect(ontology, go_term):
    if not go_term.startswith("GO:"):
        return None
    else:
        # Check ancestors for root terms
        if is_molecular_function(ontology, go_term):
            return 'F'
        elif is_cellular_component(ontology, go_term):
            return 'C'
        elif is_biological_process(ontology, go_term):
            return 'P'
###################################################

# Stealing a lot of code for this from ontobio.rdfgen:
# https://github.com/biolink/ontobio


def expand_uri_wrapper(id):
    uri = expand_uri(id, cmaps=[prefix_context])
    return uri

HAS_SUPPORTING_REFERENCE = URIRef(expand_uri_wrapper("dc:source"))
ENABLED_BY = URIRef(expand_uri_wrapper(ro.enabled_by))
ENABLES = URIRef(expand_uri_wrapper(ro.enables))
INVOLVED_IN = URIRef(expand_uri_wrapper(ro.involved_in))
PART_OF = URIRef(expand_uri_wrapper(ro.part_of))
OCCURS_IN = URIRef(expand_uri_wrapper(ro.occurs_in))
COLOCALIZES_WITH = URIRef(expand_uri_wrapper(ro.colocalizes_with))
CONTRIBUTES_TO = URIRef(expand_uri_wrapper("RO:0002326"))
MOLECULAR_FUNCTION = URIRef(expand_uri_wrapper(upt.molecular_function))
REGULATES = URIRef(expand_uri_wrapper("RO:0002211"))

now = datetime.datetime.now()


class Annoton():
    def __init__(self, subject_id, assocs, connections=None):
        self.enabled_by = subject_id
        self.annotations = assocs
        self.connections = connections
        self.individuals = {}

class GoCamModel():
    relations_dict = {
        "has_direct_input": "RO:0002400",
        "has input": "RO:0002233",
        "has_regulation_target": "RO:0002211",  # regulates
        "regulates_activity_of": "RO:0002578",  # directly regulates
        "with_support_from": "RO:0002233",  # has input
        "directly_regulates": "RO:0002578",
        "directly_positively_regulates": "RO:0002629",
        "directly_negatively_regulates": "RO:0002630"
    }

    def __init__(self, modeltitle, connection_relations=None):
        cam_writer = CamTurtleRdfWriter(modeltitle)
        self.writer = AnnotonCamRdfTransform(cam_writer)
        self.modeltitle = modeltitle
        self.classes = []
        self.individuals = {}   # Maintain entity-to-IRI dictionary. Prevents dup individuals but we may want dups?
        self.graph = networkx.MultiDiGraph()  # networkx graph of individuals and relations? Could this replace self.individuals? Will this conflict with self.writer.writer.graph?
        # Each node:
        ## node_id
        ## class
        ## attributes
        # Each edge:
        ## source
        ## target
        ## relation
        ## other attributes?
        if connection_relations is None:
            self.connection_relations = GoCamModel.relations_dict
        else:
            self.connection_relations = connection_relations
        self.declare_properties()

    def write(self, filename):
        if path.splitext(filename)[1] != ".ttl":
            filename += ".ttl"
        with open(filename, 'wb') as f:
            self.writer.writer.serialize(destination=f)

    def declare_properties(self):
        # AnnotionProperty
        self.writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)

    def declare_class(self, class_id):
        if class_id not in self.classes:
            self.writer.emit_type(URIRef("http://identifiers.org/" + class_id), OWL.Class)
            self.classes.append(class_id)

    def declare_individual(self, entity_id):
        entity = genid(base=self.writer.writer.base + '/')
        self.writer.emit_type(entity, self.writer.uri(entity_id))
        self.writer.emit_type(entity, OWL.NamedIndividual)
        self.individuals[entity_id] = entity
        self.graph.add_node(entity, **{"label": entity_id})
        return entity

    def add_axiom(self, statement, evidence=None):
        (source_id, property_id, target_id) = statement
        stmt_id = self.find_bnode(statement)
        if stmt_id is None:
            stmt_id = self.writer.blanknode()
            self.writer.emit_type(stmt_id, OWL.Axiom)
        self.writer.emit(stmt_id, OWL.annotatedSource, source_id)
        self.writer.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.writer.emit(stmt_id, OWL.annotatedTarget, target_id)

        if evidence:
            self.add_evidence(stmt_id, evidence.evidence_code, evidence.references)

        return stmt_id

    def add_evidence(self, axiom, evidence_code, references):
        ev = GoCamEvidence(evidence_code, references)
        # Try finding existing evidence object containing same type and references
        # ev_id = self.writer.find_or_create_evidence_id(ev)
        ev_id = self.writer.create_evidence(ev)
        self.writer.emit(axiom, URIRef("http://geneontology.org/lego/evidence"), ev_id)

    def add_connection(self, gene_connection, source_annoton):
        # Switching from reusing existing activity node from annoton to creating new one for each connection - Maybe SPARQL first to check if annoton activity already used for connection?
        # Check annoton for existing activity.
        # if gene_connection.object_id in source_annoton.individuals:
        #     # If exists and activity has connection relation,
        #     # Look for two triples: (gene_connection.object_id, ENABLED_BY, source_annoton.enabled_by) and (gene_connection.object_id, connection_relations, anything)
        # Annot MF should be declared by now - don't declare object_id if object_id == annot MF?
        if gene_connection.gp_b not in self.individuals:
            return
        source_id = None
        uri_list = self.uri_list_for_individual(gene_connection.object_id)
        for u in uri_list:
            if gene_connection.relation in self.connection_relations:
                rel = URIRef(expand_uri_wrapper(self.connection_relations[gene_connection.relation]))
                # Annot MF should be declared by now - don't declare object_id if object_id == annot MF?
                try:
                    annot_mf = source_annoton.molecular_function["object"]["id"]
                except:
                    annot_mf = ""
                if self.writer.writer.graph.__contains__((u,rel,None)) and gene_connection.object_id != annot_mf:
                    source_id = self.declare_individual(gene_connection.object_id)
                    source_annoton.individuals[gene_connection.object_id] = source_id
                    break

        if source_id is None:
            try:
                source_id = source_annoton.individuals[gene_connection.object_id]
            except KeyError:
                source_id = self.declare_individual(gene_connection.object_id)
                source_annoton.individuals[gene_connection.object_id] = source_id
        # Add enabled by stmt for object_id - this is essentially adding another annoton connecting gene-to-extension/with-MF to the model
        self.writer.emit(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        self.writer.emit_axiom(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        property_id = URIRef(expand_uri_wrapper(self.connection_relations[gene_connection.relation]))
        target_id = self.individuals[gene_connection.gp_b]
        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
        self.writer.emit(source_id, property_id, target_id)
        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
        self.writer.emit_axiom(source_id, property_id, target_id)

    def uri_list_for_individual(self, individual):
        uri_list = []
        graph = self.writer.writer.graph
        for t in graph.triples((None,None,self.writer.uri(individual))):
            uri_list.append(t[0])
        return uri_list

    def triples_by_ids(self, subject, relation_uri, object_id):
        graph = self.writer.writer.graph

        triples = []
        if subject.__class__.__name__ == "URIRef" or subject is None:
            subjects = [subject]
        else:
            subjects = self.uri_list_for_individual(subject)
        if object_id.__class__.__name__ == "URIRef" or object_id is None:
            objects = [object_id]
        else:
            objects = self.uri_list_for_individual(object_id)
        for object_uri in objects:
            for subject_uri in subjects:
                # if (subject_uri, relation_uri, object_uri) in graph:
                #     triples.append((subject_uri, relation_uri, object_uri))
                for t in graph.triples((subject_uri, relation_uri, object_uri)):
                    triples.append(t)
        return triples

    def individual_label_for_uri(self, uri):
        ind_list = []
        graph = self.writer.writer.graph
        for t in graph.triples((uri, RDF.type, None)):
            if t[2] != OWL.NamedIndividual: # We know OWL.NamedIndividual triple does't contain the label so don't return it
                ind_list.append(t[2])
        return ind_list

    def axioms_for_source(self, source, property_uri=None):
        if property_uri is None:
            property_uri = OWL.annotatedSource
        axiom_list = []
        graph = self.writer.writer.graph
        for uri in self.uri_list_for_individual(source):
            for t in graph.triples((None, property_uri, uri)):
                axiom_list.append(t[0])
        return axiom_list

    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]


class AssocGoCamModel(GoCamModel):
    def __init__(self, modeltitle, assocs, ontology, connection_relations=None):
        GoCamModel.__init__(self, modeltitle, connection_relations)
        self.associations = assocs
        self.ontology = ontology

    def translate(self):
        for a in self.associations:

            annoton = Annoton(a["subject"]["id"], [a])

            self.declare_class(annoton.enabled_by)
            gp_uri = self.declare_individual(annoton.enabled_by)
            term = a["object"]["id"]
            self.declare_class(term)
            term_uri = self.declare_individual(term)

            # Paul's current rules are based on aspect, similar to rdfgen's current state, which may change
            # since relation is explicitly stated in GPAD
            # Standardize aspect using GPAD relations?

            ### Use these commented lines instead once https://github.com/biolink/ontobio/pull/281 is in ontobio ###
            # aspector = GoAspector(self.ontology)
            # aspect = aspector.go_aspect(term)
            aspect = go_aspect(self.ontology, term)

            aspect_triples = []
            # Axiom time! - Stealing from ontobio/rdfgen
            if aspect == 'F':
                aspect_triples.append(self.writer.emit(term_uri, ENABLED_BY, gp_uri))
            elif aspect == 'P':
                self.declare_class(upt.molecular_function)
                mf_root_uri = self.declare_individual(upt.molecular_function)
                aspect_triples.append(self.writer.emit(mf_root_uri, ENABLED_BY, gp_uri))
                aspect_triples.append(self.writer.emit(mf_root_uri, PART_OF, term_uri))
            elif aspect == 'C':
                self.declare_class(upt.molecular_function)
                mf_root_uri = self.declare_individual(upt.molecular_function)
                aspect_triples.append(self.writer.emit(mf_root_uri, ENABLED_BY, gp_uri))
                aspect_triples.append(self.writer.emit(mf_root_uri, OCCURS_IN, term_uri))

            # Add evidence
            for atr in aspect_triples:
                axiom_id = self.add_axiom(atr)
                self.add_evidence(axiom_id, a["evidence"]["type"],
                               a["evidence"]["has_supporting_reference"])

            # Translate extension - maybe add function argument for custom translations?


class GoCamEvidence():
    def __init__(self, code, references):
        self.evidence_code = code
        self.references = references
        self.id = None

class CamTurtleRdfWriter(TurtleRdfWriter):
    def __init__(self, modeltitle):
        self.base = genid(base="http://model.geneontology.org")
        self.graph = rdflib.Graph(identifier=self.base)
        self.graph.bind("owl", OWL)
        self.graph.bind("obo", "http://purl.obolibrary.org/obo/")
        self.graph.bind("dc", DC)

        self.graph.add((self.base, RDF.type, OWL.Ontology))

        # Model attributes TODO: Should move outside init
        self.graph.add((self.base, URIRef("http://purl.org/pav/providedBy"), Literal("http://geneontology.org")))
        self.graph.add((self.base, DC.date, Literal(str(now.year) + "-" + str(now.month) + "-" + str(now.day))))
        self.graph.add((self.base, DC.title, Literal(modeltitle)))
        self.graph.add((self.base, DC.contributor, Literal("http://orcid.org/0000-0002-6659-0416"))) #TODO
        self.graph.add((self.base, URIRef("http://geneontology.org/lego/modelstate"), Literal("development")))
        self.graph.add((self.base, OWL.versionIRI, self.base))
        self.graph.add((self.base, OWL.imports, URIRef("http://purl.obolibrary.org/obo/go/extensions/go-lego.owl")))

class AnnotonCamRdfTransform(CamRdfTransform):
    def __init__(self, writer=None):
        CamRdfTransform.__init__(self, writer)
        self.annotons = []
        self.classes = []
        self.evidences = []
        self.ev_ids = []
        self.bp_id = None

    # TODO Remove "find" feature
    def find_or_create_evidence_id(self, evidence):
        for existing_evidence in self.evidences:
            if evidence.evidence_code == existing_evidence.evidence_code and set(evidence.references) == set(existing_evidence.references):
                if existing_evidence.id is None:
                    existing_evidence.id = genid(base=self.writer.base + '/')
                    self.ev_ids.append(existing_evidence.id)
                return existing_evidence.id
        return self.create_evidence(evidence)

    def create_evidence(self, evidence):
        # Use/figure out standard for creating URIs
        # Find minerva code to generate URI, add to Noctua doc
        ev_id = genid(base=self.writer.base + '/')
        evidence.id = ev_id
        # ev_cls = self.eco_class(self.uri(evidence.evidence_code))
        # ev_cls = self.eco_class(evidence.evidence_code) # This is already ECO:##### due to a GPAD being used
        ev_cls = self.uri(evidence.evidence_code)
        self.emit_type(ev_id, OWL.NamedIndividual)
        self.emit_type(ev_id, ev_cls)
        for ref in evidence.references:
            o = Literal(ref) # Needs to go into Noctua like 'PMID:####' rather than full URL
            self.emit(ev_id, HAS_SUPPORTING_REFERENCE, o)
        self.evidences.append(evidence)
        return evidence.id

    # Use only for OWLAxioms
    # There are two of these methods. AnnotonCamRdfTransform.find_bnode and GoCamModel.find_bnode. Which one is used?
    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]

    def emit_axiom(self, source_id, property_id, target_id):
        stmt_id = self.blanknode()
        self.emit_type(stmt_id, OWL.Axiom)
        self.emit(stmt_id, OWL.annotatedSource, source_id)
        self.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.emit(stmt_id, OWL.annotatedTarget, target_id)
        return stmt_id

    def find_annotons(self, enabled_by, annotons_list=None):
        found_annotons = []
        if annotons_list is not None:
            annotons = annotons_list
        else:
            annotons = self.annotons
        for annoton in annotons:
            if annoton.enabled_by == enabled_by:
                found_annotons.append(annoton)
        return found_annotons

    def add_individual(self, individual_id, annoton):
        obj_uri = self.uri(individual_id)
        if individual_id not in annoton.individuals:
            tgt_id = genid(base=self.writer.base + '/')
            annoton.individuals[individual_id] = tgt_id
            self.emit_type(tgt_id, obj_uri)
            self.emit_type(tgt_id, OWL.NamedIndividual)
        else:
            tgt_id = annoton.individuals[individual_id]
