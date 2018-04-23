"""
Get aggregate counts of node labels per edge type
Note that nodes contain a set of labels; therefore,
we are counting the number of nodes with a certain
set of labels eg
versionIRI	outgoing	17212	['sequence feature', 'cliqueLeader', 'Node']
versionIRI	outgoing	4546	['cliqueLeader', 'Node']
versionIRI	incoming	280488	['Class', 'gene', 'sequence feature', 'cliqueLeader', 'Node']
versionIRI	incoming	6135	['cliqueLeader', 'Class', 'gene', 'molecular entity', 'sequence feature', 'Node']
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
    MATCH (n)-[:`http://purl.obolibrary.org/obo/GENO_0000418`]->()
    RETURN count(labels(n)) as count, labels(n) as labels
    ORDER BY count DESC
    LIMIT 5
    """.format(rel)
    for aggregate in get_scigraph_results(label_query):
        print("{}\toutgoing\t{}\t{}".format(rel, aggregate['count'], aggregate['labels']))
    label_query = """
        MATCH (n)<-[:`http://purl.obolibrary.org/obo/GENO_0000418`]-()
        RETURN count(labels(n)) as count, labels(n) as labels
        ORDER BY count DESC
        LIMIT 5
        """.format(rel)
    for aggregate in get_scigraph_results(label_query):
        print("{}\tincoming\t{}\t{}".format(rel, aggregate['count'], aggregate['labels']))
