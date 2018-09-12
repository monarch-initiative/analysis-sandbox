"""
Get aggregate counts of node labels per edge type
Note that nodes contain a set of labels; therefore,
we are counting the number of nodes with a certain
set of labels eg
isDefinedBy	outgoing	2112891	['Class', 'gene', 'sequence feature', 'cliqueLeader', 'Node']
isDefinedBy	outgoing	1024789	['NamedIndividual', 'publication', 'cliqueLeader', 'Node']
isDefinedBy	incoming	5613517	['Ontology', 'cliqueLeader', 'Node']
isDefinedBy	incoming	353	['Ontology', 'NamedIndividual', 'cliqueLeader', 'Node']
"""

from neo4j.v1 import GraphDatabase

query = "call db.relationshipTypes"
rel_types = []

scigraph_bolt = "bolt://neo4j.monarchinitiative.org:443"
driver = GraphDatabase.driver(scigraph_bolt, auth=("neo4j", "neo4j"))

def get_scigraph_results(query):
    with driver.session() as session:
        with session.begin_transaction() as tx:
            for record in tx.run(query):
                yield record

for record in get_scigraph_results(query):
    rel_types.append(record['relationshipType'])

for rel in rel_types:
    label_query = """
    MATCH (n)-[:`{}`]->()
    RETURN count(labels(n)) as count, labels(n) as labels
    ORDER BY count DESC
    LIMIT 5
    """.format(rel)
    for aggregate in get_scigraph_results(label_query):
        print("{}\toutgoing\t{}\t{}".format(rel, aggregate['count'], aggregate['labels']))
    label_query = """
        MATCH (n)<-[:`{}`]-()
        RETURN count(labels(n)) as count, labels(n) as labels
        ORDER BY count DESC
        LIMIT 5
        """.format(rel)
    for aggregate in get_scigraph_results(label_query):
        print("{}\tincoming\t{}\t{}".format(rel, aggregate['count'], aggregate['labels']))
