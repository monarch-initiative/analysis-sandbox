from ontobio.ontol_factory import OntologyFactory
from ontobio.io.ontol_renderers import OboFormatGraphRenderer
from ontobio.io.ontol_renderers import GraphRenderer
import requests
import ontobio.sparql.wikidata as wd
from ontobio.assoc_factory import AssociationSetFactory
from dipper.graph.RDFGraph import RDFGraph
from scigraph.api.SciGraph import SciGraph
import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Make ontology from wikidata

onto_factory = OntologyFactory()
wd_ontology = onto_factory.create('wdq:Q185034')  # Sickle cell anemia
qids = wd_ontology.search('Sickle%')

# Traverse up and down from query node in our sub-ontology
nodes = wd_ontology.traverse_nodes(qids, up=True, down=True)

renderer = GraphRenderer.create('obo')
renderer.outfile = './output/wd-ontology.obo'
# renderer.write(wd_ontology)
# >> AttributeError: 'EagerWikidataOntology' object has no attribute 'all_logical_definitions'

renderer.write_subgraph(wd_ontology, nodes, query_ids=qids)

# Get GO terms
outfile = open('./output/go-terms.tsv', 'w')

[ptsd] = wd_ontology.search('Sickle Cell Anemia')
proteins = wd.canned_query('disease2protein', ptsd)

go = onto_factory.create('go')

afactory = AssociationSetFactory()
aset = afactory.create(ontology=go,
                       subject_category='gene',
                       object_category='function',
                       taxon='NCBITaxon:9606')

for n in wd_ontology.nodes():
    proteins = wd.canned_query('disease2protein', n)
    anns = [a for p in proteins for a in aset.annotations(p)]
    if len(anns) > 0:
        print("{} {}".format(n, wd_ontology.label(n)))
        for a in anns:
            outfile.write("{}\t{}\n".format(a, go.label(a)))

# Endpoints
SCIGRAPH_ONTOLOGY = 'https://scigraph-ontology-dev.monarchinitiative.org/scigraph/'
SCIGRAPH_DATA = 'https://scigraph-data-dev.monarchinitiative.org/scigraph/'
GOLR_URL = 'https://solr.monarchinitiative.org/solr/golr/select'

# Get mondo subset
sickle_cell_anemia = 'OMIM:603903'
disease_graph = RDFGraph()

sg = SciGraph(SCIGRAPH_ONTOLOGY)
parent_graph = sg.neighbors(sickle_cell_anemia,
    {'depth':25, 'direction': 'OUTGOING', 'relationshipType': 'subClassOf'})
child_graph = sg.neighbors(sickle_cell_anemia,
    {'depth':2, 'direction': 'INCOMING', 'relationshipType': 'subClassOf'})
eq_graph = sg.neighbors(sickle_cell_anemia, {'depth':10, 'relationshipType': 'equivalentClass'})


def add_triples_from_bbop(bbop_graph, rdf_graph):
    for e in bbop_graph.edges:
        if not re.search(r':', e.predicate):
            if e.predicate == 'subClassOf':
                e.predicate = 'rdfs:subClassOf'
            elif e.predicate == 'equivalentClass':
                e.predicate = 'owl:equivalentClass'
            elif e.predicate == 'isDefinedBy':
                e.predicate = 'rdfs:isDefinedBy'

        rdf_graph.addTriple(e.subject, e.predicate, e.target)
    return rdf_graph

disease_graph = add_triples_from_bbop(parent_graph, disease_graph)
disease_graph = add_triples_from_bbop(child_graph, disease_graph)
disease_graph = add_triples_from_bbop(eq_graph, disease_graph)

mondo_stub = './output/mondo-stub.xml'
mondo_stub_ttl = './output/mondo-stub.ttl'

disease_graph.serialize(mondo_stub, 'xml')
disease_graph.serialize(mondo_stub_ttl, 'ttl')

data_graph = RDFGraph()
sg = SciGraph(SCIGRAPH_DATA)


# Get all children + sickle cell
children = [node.id for node in child_graph.nodes if node.id != "owl:Nothing"]

for node in children:
    has_phenotype = sg.neighbors(node,
        {'depth':1, 'direction': 'OUTGOING', 'relationshipType': 'RO:0002200'})
    data_graph = add_triples_from_bbop(has_phenotype, data_graph)

data_graph.addTriple('NCBIGene:3043', 'RO:0002200', sickle_cell_anemia)

filters = [
    ['object_closure: "OMIM:603903"', 'subject_category:"variant"', 'object_category:"disease"'],
    ['object_closure: "OMIM:603903"','subject_category:"model"', 'object_category:"disease"'],
    ['subject_closure: "OMIM:603903"','subject_category:"disease"', 'object_category:"pathway"'],
    ['subject_closure: "OMIM:603903"','subject_category:"disease"', 'object_category:"phenotype"'],
]

for filtr in filters:
    golr_params = {
        'wt': 'json',
        'q': '*:*',
        'fq': filtr,
        'fl': 'subject, subject_label, subject_gene, relation, object, object_label',
        'rows': 1000,
        'start': 0
    }

    resultCount = golr_params['rows']

    while golr_params['start'] < resultCount:
        solr_request = requests.get(GOLR_URL, params=golr_params)
        response = solr_request.json()
        resultCount = response['response']['numFound']

        for doc in response['response']['docs']:
            if 'subject_gene' in doc:
                data_graph.addTriple(doc['subject'], 'GENO:0000418', doc['subject_gene'])
            if 'relation' in doc:
                data_graph.addTriple(doc['subject'], doc['relation'], doc['object'])
            if 'subject_label' in doc:
                data_graph.addTriple(doc['subject'], 'rdfs:label', doc['subject_label'], True)
            if 'object_label' in doc:
                data_graph.addTriple(doc['object'], 'rdfs:label', doc['object_label'], True)

        if golr_params['start'] % 100000 == 0:
            logger.info("Processed {} documents".format(golr_params['start']))

        golr_params['start'] += golr_params['rows']


data_xml = './output/data-graph.xml'
data_ttl = './output/data-graph.ttl'

data_graph.serialize(data_xml, 'xml')
data_graph.serialize(data_ttl, 'ttl')

