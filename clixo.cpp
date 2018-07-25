#include <iostream>
#include <stdlib.h>
#include "dagConstruct.h"


int main(int argc, char* argv[]) {

  if (argc < 3) {
    cout << "Needs 3 arguments: 1) file with pairwise distances between elements.  ";
    cout << "File should be 3 tab separated columns with 1 undirected edge per line of the format node1, node2, edgeWeight (similarity between node1 and node2)" << endl;
    cout << "2) threshold between clusters (alpha parameter)" << endl;
    cout << "3) merge density for overlapping clusters (beta parameter. Optional with default = 0.5)" << endl;
    cout << "4) parameter regarding Newman's modularity. Optional with default = 0.002)" << endl;
    cout << "5) if given entering debug mode, print more information" << endl;
    return 0;
  }

  /*for (unsigned i = 0; i < 600; ++i) {
    cout << i << "\t" << Corrector::correction(i) << endl;
  }
  return 1;*/

  map<string, unsigned> nodeNamesToIDs;

  // Load input network graphs
  string netFile = argv[1];
  double threshold = stod(argv[2]);
  double density = 0.5;
  double modular = 0.25;
  bool debug = false;

  if (argc >= 4) {
    density = stod(argv[3]);
  }

  if (argc >= 5) {
    modular = stod(argv[4]);
  }

  if (argc >= 6) {
    debug = true;
  }
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
  dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances, threshold, density, modular, debug);
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
