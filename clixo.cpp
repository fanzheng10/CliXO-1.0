#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "dagConstruct.h"

void usage(char *prog_name) {
    cout << "USAGE: " << prog_name << " -i input_file -a alpha [-b beta] [-m Newman's modularity cutoff] [-z z-modularity cutoff] [-s stop_score]" << endl << endl;
    cout << "-i\tfile with pairwise distances between elements. ";
    cout << "File should be 3 tab separated columns with 1 undirected edge per line of the format node1, node2, edgeWeight (similarity between node1 and node2)" << endl;
    cout << "-a\tthreshold between clusters (alpha parameter)" << endl;
    cout << "-b\t(optional, default = 0.5) merge density for overlapping clusters (beta parameter)" << endl;
    cout << "-m\t(optional, default = 0.005) cutoff for filtering clusters based on Newman's modularity" << endl;
    cout << "-z\t(optional, default = 0.2) cutoff for filtering clusters based on z-score modularity" << endl;
    cout << "-s\t(optional, default = 0) score threshold to stop the program" << endl;
    exit(0);
}

int main(int argc, char* argv[]) {

    //default values
    string netFile;
    bool netflag=false;
    double alpha = 0;
    bool aflag=false;
    double beta = 0.5;
    double modular = 0.2;
    double zmodular = 0.005;
    double stopt = 0;

    //parse argument
    int c;
    while ((c = getopt(argc, argv, "i:a:b:m:z:s:")) != -1) {
        switch (c) {
            case 'i':
                netFile = optarg;
                netflag=true;
                break;
            case 'a':
                alpha = stod(optarg);
                aflag=true;
                break;
            case 'b':
                beta = stod(optarg);
                break;
            case 'm':
                modular = stod(optarg);
                break;
            case 'z':
                zmodular = stod(optarg);
                break;
            case 's':
                stopt = stod(optarg);
                break;
            case '?':
                fprintf (stderr, "Unknown option `-%c' \n", optopt);
                usage(argv[0]);
                break;
        }
    }
    if (!netflag) {
        fprintf(stderr, "-i is mandatory\n");
        usage(argv[0]);
    }
    if (!aflag) {
        fprintf(stderr, "-a is mandatory\n");
        usage(argv[0]);
    }

    map<string, unsigned> nodeNamesToIDs;

    string terminalName = "gene";
    time_t start, end;
    time(&start);
    cout << "# Loading input network graph" << endl;
    graph_undirected inputNetwork(netFile, nodeNamesToIDs);
    time(&end);
    double dif = difftime(end,start);
    cout << "# Loading input network took " << dif << " seconds" << endl;

    // Create ontology from other networks
    cout << "# Clique finding beginning" << endl;
    DAGraph ontology;
    ontology.setTerminalName(terminalName);
    nodeDistanceObject nodeDistances;

    time (&start);
    dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances, alpha, beta, modular, zmodular, stopt);
    time (&end);
    dif = difftime(end,start);
    cout << "# Ontology construction took " << dif << " seconds" << endl;
    cout << "# Ontology is: " << endl;

    for(map< pair<unsigned,unsigned>, string >::iterator edgesIt = ontology.edgesBegin(); edgesIt != ontology.edgesEnd();++edgesIt) {
        // THIS VERSION GIVES THE DISTANCE BETWEEN POINTS
        //cout << ontology.getName(edgesIt->first.first) << "\t" << ontology.getName(edgesIt->first.second) << "\t" << edgesIt->second << "\t" << (ontology.getWeight(edgesIt->first.second) - ontology.getWeight(edgesIt->first.first)) / 2.0 << endl;

        // THIS VERSION JUST TELLS THE WEIGHT ON THE PARENT TERM
        cout << ontology.getName(edgesIt->first.first) << "\t" << ontology.getName(edgesIt->first.second) << "\t" << edgesIt->second << "\t" << ontology.getWeight(edgesIt->first.first) << endl;
    }
    return 1;
}
