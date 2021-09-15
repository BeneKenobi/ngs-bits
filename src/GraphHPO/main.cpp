#include "ToolBase.h"
#include "Graph.h"
#include "OntologyTermCollection.h"
#include "Exceptions.h"
#include "Helper.h"
#include <QByteArray>
#include <QString>
#include <QHash>
#include <QSet>
#include <QList>
#include <QStringList>
#include <QTextStream>
#include <cmath>


struct NodeContent
{
    int associated_diseases;
    int directly_associated_diseases;
    double information_content;

    NodeContent()
        : associated_diseases(0),
          directly_associated_diseases(0),
          information_content(0.0)
    {
    }
};

struct EdgeContent
{
    double weight;

    EdgeContent()
        : weight(0.0)
    {
    }
};

class ConcreteTool
        : public ToolBase
{
    Q_OBJECT

private:
    QHash<QString, QStringList> annotations;
    QHash<QString, QStringList> associated_terms;
    Graph<NodeContent, EdgeContent> hpo_graph;

    // create a directed graph from the ontology term collection
    Graph<NodeContent, EdgeContent> parseTermCollection(OntologyTermCollection &terms)
    {
        Graph<NodeContent, EdgeContent> hpo_graph(true);

        for(int i = 0; i < terms.size(); i++)
        {
            OntologyTerm term = terms.get(i);

            foreach(QByteArray parent, term.parentIDs())
            {
                NodeContent node_content_1;
                NodeContent node_content_2;
                EdgeContent edge_content;
                hpo_graph.addEdge(parent, node_content_1, term.id(), node_content_2, edge_content);
            }
        }

        return hpo_graph;
    }

    // parse diseases associated with each HPO term
    QHash<QString, QStringList> parseAnnotations(const QString &file)
    {
        QHash<QString, QStringList> annotations;

        QSharedPointer<QFile> reader = Helper::openFileForReading(file);
        QTextStream in(reader.data());

        while(!in.atEnd())
        {
            QString line = in.readLine();

            // skip comments
            if(line.startsWith("#"))
            {
                continue;
            }

            // file format: 1st column: DatabaseID; 4th column: HPO_ID
            QStringList anno = line.split("\t");

            if(!annotations.contains(anno.at(3)))
            {
                annotations.insert(anno.at(3), QStringList());
            }

            // add the disease (ID) to the list of associated diseases for the HPO term
            if(!annotations[anno.at(3)].contains(anno.at(0)))
            {
                annotations[anno.at(3)].append(anno.at(0));
            }
        }

        return annotations;
    }

    // parse HPO terms associated with each gene
    QHash<QString, QStringList> parseAssociatedPhenotypes(const QString &file)
    {
        QHash<QString, QStringList> associated_terms;

        QSharedPointer<QFile> reader = Helper::openFileForReading(file);
        QTextStream in(reader.data());

        while(!in.atEnd())
        {
            QString line = in.readLine();

            // skip comments
            if(line.startsWith("#"))
            {
                continue;
            }

            // file format: 1st column: entrez-gene-id; 3rd column: HPO-Term-ID
            QStringList entry = line.split("\t");

            if(!associated_terms.contains(entry.at(0)))
            {
                associated_terms.insert(entry.at(0), QStringList());
            }

            if(!associated_terms[entry.at(0)].contains(entry.at(2)))
            {
                associated_terms[entry.at(0)].append(entry.at(2));
            }
        }

        return associated_terms;
    }

    // recursively count number of directly and indirectly associated diseases for a term
    QSet<QString> countAssociatedDiseases(QString node_name)
    {
        QSet<QString> associated_diseases = QSet<QString>::fromList(annotations[node_name]);
        hpo_graph.getNode(node_name).data()->nodeContent().directly_associated_diseases = associated_diseases.size();

        Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, hpo_graph.getAdjacentNodes(node_name))
        {
            associated_diseases.unite(countAssociatedDiseases(node.data()->nodeName()));
        }

        hpo_graph.getNode(node_name).data()->nodeContent().associated_diseases = associated_diseases.size();

        return associated_diseases;
    }

    // find the most informative common ancestor (MICA) of two terms and return its information content
    double findMICA(QString term_1, QString term_2, OntologyTermCollection& terms)
    {
        if(!hpo_graph.hasNode(term_1) || !hpo_graph.hasNode(term_2))
        {
            THROW(ArgumentException, "At least one term is not contained in the graph.");
        }

        QSet<QString> parents_term_1;
        parents_term_1.insert(term_1);
        QSet<QString> parents_term_2;
        parents_term_2.insert(term_2);
        QSet<QString> current_terms_1;
        current_terms_1.insert(term_1);
        QSet<QString> current_terms_2;
        current_terms_2.insert(term_2);

        //QTextStream out(stdout);

        // go up the DAG until at least one common ancestor was found
        while(!parents_term_1.intersects(parents_term_2))
        {
            if(!terms.containsByID(term_1.toUtf8()) || !terms.containsByID(term_2.toUtf8()))
            {
                THROW(ArgumentException, "At least one term is not contained in the ontology");
            }

            QSet<QString> new_terms;
            foreach(QString term, current_terms_1)
            {
                foreach(QString parent, terms.getByID(term.toUtf8()).parentIDs())
                {
                    parents_term_1.insert(parent);
                    new_terms.insert(parent);
                }
            }
            current_terms_1 = new_terms;

            new_terms.clear();
            foreach(QString term, current_terms_2)
            {
                foreach(QString parent, terms.getByID(term.toUtf8()).parentIDs())
                {
                    parents_term_2.insert(parent);
                    new_terms.insert(parent);
                }
            }
            current_terms_2 = new_terms;
        }

        // determine the maximum information content (in case the intersection contains more than one ancestor)
        double max_information_content = 0.0;
        foreach(QString term, parents_term_1.intersect(parents_term_2))
        {
            double information_content = hpo_graph.getNode(term).data()->nodeContent().information_content;
            //out << term_1 << "\t" << term_2 << "\t" << term << "\t" << information_content << endl;
            if(information_content > max_information_content)
            {
                max_information_content = information_content;
            }
        }

        return max_information_content;
    }

    // determine the score for a gene given a list of HPO-terms (maximum similarity)
    double scoreGene(QString gene_id, QStringList hpo_terms, OntologyTermCollection& terms)
    {
        if(!associated_terms.contains(gene_id))
        {
            return 0.0;
        }

        QList<double> scores;
        //QTextStream out(stdout);

        foreach(QString gene_term, associated_terms[gene_id])
        {
            // terms associated with genes could also be from another category, not just from phenotypic abnormality
            if(!hpo_graph.hasNode(gene_term))
            {
                continue;
            }

            double max_score = 0.0;

            // take the phenotype term with maximum similarity with the current gene-associated term
            foreach(QString hpo_term, hpo_terms)
            {
                double score = findMICA(gene_term, hpo_term, terms);
                if(score > max_score)
                {
                    max_score = score;
                }
            }
            scores.append(max_score);
            //out << gene_term << "\t" << max_score << endl;
        }

        // average over the scores for all terms associated with the gene
        double sum = 0.0;
        foreach(double score, scores)
        {
            sum += score;
        }

        return sum / scores.size();
    }

public:
    ConcreteTool(int& argc, char *argv[])
        : ToolBase(argc, argv)
    {
    }

    virtual void setup()
    {
        setDescription("Creates simple representation of HPO directed acyclic graph.");
        addInfile("in", "Input TSV file with comma-separated list of HPO terms and gene (Entrez Gene ID)", false);
        addInfile("obo", "Human phenotype ontology obo file.", false);
        addInfile("annotation", "TSV file with HPO annotations to diseases", false);
        addInfile("associated-terms", "TSV file with all associated HPO terms for each gene", false);
        addOutfile("out", "Output TSV file with a score for each HPO-term list/gene combination", false);
    }

    virtual void main()
    {
        // init
        OntologyTermCollection terms(getInfile("obo"), true);
        hpo_graph = parseTermCollection(terms).extractSubgraph("HP:0000118");

        annotations = parseAnnotations(getInfile("annotation"));

        associated_terms = parseAssociatedPhenotypes(getInfile("associated-terms"));

        // recursively count number of associated diseases, starting from root node "Phenotypic abnormality"
        int total_diseases = countAssociatedDiseases("HP:0000118").size();

        // determine information content of each node
        Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, hpo_graph.adjacencyList().keys())
        {
            if(node.data()->nodeContent().associated_diseases != 0)
            {
                node.data()->nodeContent().information_content =
                        -log((double) node.data()->nodeContent().associated_diseases / total_diseases);
            }
            else
            {
                node.data()->nodeContent().information_content = 0;
            }
        }

        //hpo_graph.store(getOutfile("out"));

        /*QTextStream out(stdout);
        out << findMICA("HP:0000494", "HP:0000506", terms) << endl;*/

        QSharedPointer<QFile> writer = Helper::openFileForWriting(getOutfile("out"));
        QTextStream stream(writer.data());
        stream << "#hpo_terms\tgene_id\tscore" << endl;

        QSharedPointer<QFile> reader = Helper::openFileForReading(getInfile("in"));
        QTextStream in(reader.data());

        while(!in.atEnd())
        {
            QString line = in.readLine();

            // skip comments
            if(line.startsWith("#"))
            {
                continue;
            }

            QStringList entry = line.split("\t");
            double score = scoreGene(entry.at(1), entry.at(0).split(","), terms);
            stream << line << "\t" << score << endl;
        }
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
    ConcreteTool tool(argc, argv);
    return tool.execute();
}

