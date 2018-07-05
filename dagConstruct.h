//
// Created by fanzheng on 3/20/18.
//

#ifndef DAGCONSTRUCT
#define DAGCONSTRUCT

#include <math.h>
#include <time.h>
#include <list>
#include <algorithm>
#include <thread>
#include <mutex>
#include "dag.h"
#include "graph_undirected.h"
#include "graph_undirected_bitset.h"
#include "util.h"
#include "nodeDistanceObject.h"
#include "boost/dynamic_bitset/dynamic_bitset.hpp"

bool debug = false;

// to print the gene names in the cluster
void printCluster(const boost::dynamic_bitset<unsigned long> & cluster, vector<string> & nodeIDsToNames) {
    for (unsigned i = 0; i < cluster.size(); ++i) {
        if (cluster[i]) {
            cout << nodeIDsToNames[i] << ",";
        }
    }
}

bool compPairSecondAscending(const pair< unsigned, unsigned > & i,const pair< unsigned, unsigned > & j) {
    return (i.second < j.second);
}

bool compPairSecondDescending(const pair< unsigned, unsigned > & i,const pair< unsigned, unsigned > & j) {
    return (i.second > j.second);
}

bool compClustersToCombine(const pair<pair<unsigned long, unsigned long>, double> & i, const pair<pair<unsigned long, unsigned long>, double> & j) {
    if (i.second == j.second) {
        if (i.first.first == j.first.first) {
            return i.first.second < j.first.second;
        }
        return i.first.first < j.first.first;
    }
    return i.second > j.second;
}

double calculateClusterWeight(boost::dynamic_bitset<unsigned long> cluster, nodeDistanceObject & nodeDistances) {
//        vector <double> edgeWeights;
    double sumEdgeWeights = 0;
    unsigned nEdges = cluster.count() * (cluster.count()-1)/2;
    for (unsigned long i = cluster.find_first(); i < cluster.size()-1; i=cluster.find_next(i) ) {
        for (unsigned long j = cluster.find_next(i); j < cluster.size(); j = cluster.find_next(j)) {
//                edgeWeights.push_back(nodeDistances.getDistance(i, j));
            sumEdgeWeights = sumEdgeWeights + nodeDistances.getDistance(i, j);
        }
    }
    double edgeWeightsMean = sumEdgeWeights / nEdges;
    return edgeWeightsMean;
}

class validClusterBitset {
public:
    validClusterBitset(const boost::dynamic_bitset<unsigned long> & cluster, unsigned clustID, double thisWeight) {
        elements = cluster;
        ID = clustID;
        weight = thisWeight;
        numElementsHere = cluster.count();
    }

    validClusterBitset(const vector<unsigned> & cluster, unsigned clustID, double thisWeight, unsigned numNodes) {
        elements = boost::dynamic_bitset<unsigned long>(numNodes);
        for (vector<unsigned>::const_iterator it = cluster.begin(); it != cluster.end(); ++it) {
            elements[*it] = 1;
        }
        ID = clustID;
        weight = thisWeight;
        numElementsHere = elements.count();
    }

    bool operator<(const validClusterBitset & b) const {
        return numElementsHere < b.numElements();
    }

    inline unsigned getID() {
        return ID;
    }

    inline void setID(unsigned newID) {
        ID = newID;
        return;
    }

    inline double getWeight() {
        return weight;
    }

    inline const boost::dynamic_bitset<unsigned long> & getElements() {
        return elements;
    }

    inline unsigned isElement(unsigned i) {
        return elements[i];
    }

    inline unsigned numElements() const {
        return numElementsHere;
    }

    inline void addElement(unsigned newElem) {
        elements[newElem] = 1;
        ++numElementsHere;
        return;
    }

    inline vector<unsigned> getElementsVector() {
        vector<unsigned> result;
        result.reserve(numElements());
        for (unsigned i = elements.find_first(); i < elements.size(); i = elements.find_next(i)) {
            result.push_back(i);
        }
        return result;
    }

private:
    boost::dynamic_bitset<unsigned long> elements;
    unsigned ID;
    double weight;
    unsigned numElementsHere;
};

class ClusterBitset {
public:
    ClusterBitset(const boost::dynamic_bitset<unsigned long> & cluster, unsigned long & clustID, double thisClusterWeight = 0) {
        elements = cluster;
        isClusterNew = true;
        isClusterAddedToExplain = false; // comment things that may be uncessary later
        isClusterUnexplainedCounted = false;
        ID = clustID;
        numElementsHere = cluster.count();
        clusterWeight = thisClusterWeight;
        valid = false;
        uniquelyExplainedEdges = 0;
        active = true;
    }

    ClusterBitset() {
        isClusterNew = false;
        isClusterAddedToExplain = false;
        isClusterUnexplainedCounted = false;
        ID = 0;
        numElementsHere = 0;
        clusterWeight = 0;
        valid = false;
        uniquelyExplainedEdges = 0;
        active = false;
    }

    inline double getWeight() {
        return clusterWeight;
    }

    // clusterweight minus alpha
    inline double getThresh(double alpha, unsigned numUnexplainedEdges) {
        if (numUnexplainedEdges ==0) {
            return 0;
        }
        return clusterWeight - alpha;
    }

    inline void setWeight(const double & weight) {
        clusterWeight = weight;
    }

    inline unsigned long getID() {
        return ID;
    }

    inline const boost::dynamic_bitset<unsigned long> & getElements() {
        return elements;
    }

    inline const vector<unsigned> & getElementsVector() {
        if (elementsVector.size()) {
            return elementsVector;
        }
        elementsVector.reserve(numElements());
        for (unsigned i = 0; i < elements.size(); ++i) {
            if (elements[i] == true) {
                elementsVector.push_back(i);
            }
        }
        return elementsVector;
    }

    inline unsigned numElements() const {
        return numElementsHere;
    }

    inline bool isNew() {
        return isClusterNew;
    }

    inline void setOld() {
        isClusterNew = false;
    }

    inline void setNew() {
        isClusterNew = true;
    }

    inline bool isElement(unsigned elemID) {
        return elements.test(elemID);
    }

    inline unsigned size() {
        return elements.size();
    }

    /* valid cluster is a cluster that has been printed out*/

    inline void setValid() {
        valid = true;
    }

    inline void setInvalid() {
        valid = false;
    }

    inline bool isValid() {
        return valid;
    }

    /* inactive cluster is in the middle of being considered to merge*/

    inline void setActive() {
        active = true;
    }

    inline void setInactive() {
        active = false;
    }

    inline bool isActive() {
        return active;
    }

    inline bool isAddedToExplain() {
        return isClusterAddedToExplain;
    }

    inline void setAddedToExplain() {
        isClusterAddedToExplain = true;
    }

    inline void setRemovedFromExplain() {
        isClusterAddedToExplain = false;
    }

    inline bool isUnexplainedCounted() {
        return isClusterUnexplainedCounted;
    }

    inline void setUnexplainedCounted() {
        isClusterUnexplainedCounted = true;
    }

    inline void setUnexplainedUncounted() {
        isClusterUnexplainedCounted = false;
    }

    inline unsigned getUniquelyExplainedEdges() {
        return uniquelyExplainedEdges;
    }

    inline void setUniquelyExplainedEdges(unsigned val) {
        uniquelyExplainedEdges = val;
        setUnexplainedCounted();
    }

    inline void setMergedFromID(vector <unsigned long> & ids) {
        mergedFromID = ids;
    }

    inline vector <unsigned long> getMergedFromID() {
        return mergedFromID;
    }

    inline void setMergeToID(vector <unsigned long> & ids) {
        mergedToID = ids;
    }

    inline vector <unsigned long> getMergedToID() {
        return mergedToID;
    }

    inline void addMergedFromID(unsigned long id) {
        mergedFromID.push_back(id);
    }

    inline void addMergedToID(unsigned long id) {
        mergedToID.push_back(id);
    }

    inline void removeMergedFromID(unsigned long id) {
        mergedFromID.erase(remove(mergedFromID.begin(), mergedFromID.end(), id), mergedFromID.end());
    }

    inline void removeMergedToID(unsigned long id) {
        mergedToID.erase(remove(mergedToID.begin(), mergedToID.end(), id), mergedToID.end());
    }


private:
    boost::dynamic_bitset<unsigned long> elements;
    bool isClusterNew;
    bool isClusterAddedToExplain;
    bool isClusterUnexplainedCounted;
//    bool necessary;
    unsigned long ID;
    unsigned numElementsHere;
    vector<unsigned> elementsVector;
    double clusterWeight;
    bool valid;
    unsigned uniquelyExplainedEdges;
    bool active;
    vector<unsigned long> mergedFromID; // a cluster can make other clusters non-maximal. but this cluster can be delete for some reason. When this happen, those other clusters should be reactivated.
    vector<unsigned long> mergedToID;
};



class currentClusterClassBitset {
public:
    currentClusterClassBitset(unsigned numNodesHere, double FirstWeight = 1, double thisAlpha = 0) {
        numNodes = numNodesHere;
        nodesToClusters = vector<vector<unsigned long> >(numNodes, vector<unsigned long>());
        nextID = 0;
        clustersDeleted = 0;
        clustersAdded = 0;
        newClusts = 0;
        firstWeight = FirstWeight;
        curWeight = firstWeight;
        maxNewWeight = firstWeight;
        minWeightAdded = firstWeight;
        alpha = thisAlpha;
        largestCluster = 0;

        edgesToClusters = vector<vector<unsigned> >(numNodes, vector<unsigned>(numNodes, 0));
        isEdgeExplained = vector<vector<char> >(numNodes, vector<char>(numNodes, 0));

        maxClusterSizePerNodes = vector<unsigned long>(numNodes, 0);

    };

    inline void setLargestCluster(unsigned long size) {
        largestCluster = size;
    }

    inline unsigned long getLargestCluster() {
        return largestCluster;
    }

    inline double getAlpha() {
        return alpha;
    }

    inline void setCurWeight(const double & weight) {
        if (weight < minWeightAdded) {
            minWeightAdded = weight;
        }
        curWeight = weight;
    }

    inline void setEdgeExplained(const unsigned & i, const unsigned & j) {
        isEdgeExplained[i][j] = 1;
    }

    inline bool isThisEdgeExplained(const unsigned & i, const unsigned & j) {
        return isEdgeExplained[i][j];
    }

    inline bool addClusterToExplainEdge(const unsigned & edgeEnd1, const unsigned & edgeEnd2, const unsigned long & clustID) {
        ++edgesToClusters[edgeEnd1][edgeEnd2];
        return true;
    }

    inline bool addClusterToExplanation(const vector<unsigned> & cluster, const unsigned long & clustID) {
        bool validClust = false;
        for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
            vector<unsigned>::const_iterator it2 = it1;
            ++it2;
            for ( ; it2 != cluster.end(); ++it2) {
                validClust = addClusterToExplainEdge(*it1, *it2, clustID);
            }
        }
        currentClusters[clustID].setAddedToExplain();
        return validClust;
    }

    inline void removeClusterToExplainEdge(const unsigned & edgeEnd1, const unsigned & edgeEnd2, const unsigned long & clustID) {
        --edgesToClusters[edgeEnd1][edgeEnd2];
        return;
    }

    inline void removeClusterFromExplanation(const vector<unsigned> & cluster, const unsigned long & clustID) {
        for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
            vector<unsigned>::const_iterator it2 = it1;
            ++it2;
            for ( ; it2 != cluster.end(); ++it2) {
                removeClusterToExplainEdge(*it1,*it2,clustID);
            }
        }
        currentClusters[clustID].setRemovedFromExplain();
    }

    inline void removeClusterFromExplanation(unsigned long clusterToRemove) {
        removeClusterFromExplanation(currentClusters[clusterToRemove].getElementsVector(), currentClusters[clusterToRemove].getID());
    }

    /**/
    unsigned long addCluster(const boost::dynamic_bitset<unsigned long> & newCluster,
                             nodeDistanceObject & nodeDistances,
                             vector<string> & nodeIDsToNames) {
        unsigned long newID = 0;

        // reuse some ID that are emptied due to deletion, if possible
        if (openIDs.size() != 0) {
            newID = openIDs.back();
            openIDs.pop_back();
        } else {
            newID = nextID;
            ++nextID;
        }

        if (currentClusters.size() <= newID) {
            currentClusters.resize(newID+1);
        }
        currentClusters[newID] = ClusterBitset(newCluster, newID, curWeight); // think about curWeight
        for (unsigned i = 0; i < newCluster.size() ;  ++i) {
            if (newCluster[i] == true) {
                Utils::insertInOrder(nodesToClusters[i], newID);
            }
        }
        ++clustersAdded;
        ++newClusts;
        resetClusterWeight(newID, nodeDistances);

        //check other inactivated clusters. If they are subset, add cluster relationships here
        vector <unsigned long> nonMaximal;
        nonMaximal.reserve(1000);
        for (unsigned long i = 0; i < maxClusterID(); ++i) {
            if ((numElements(i) !=0) && isNew(i) && (!isActive(i))) {
                if (getElements(i).is_subset_of(getElements(newID))) {
//                    inactivateCluster(i, nodeIDsToNames);
                    nonMaximal.push_back(i);
                }
            }
        }

        /* add relationships to the non-maximal clusters*/
        for (vector<unsigned long>::iterator nonMaximalIt = nonMaximal.begin(); nonMaximalIt < nonMaximal.end(); ++nonMaximalIt ) {
            addMergedToID(*nonMaximalIt, newID);
        }
        setMergedFromID(newID, nonMaximal);

        return newID;
    }


    inline void resetAllUnexplained() {
        for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
            if ((clustIt->numElements() != 0) && (!clustIt->isValid())) {
                clustIt->setUnexplainedUncounted();
            }
        }
    }

    inline unsigned setNumUniquelyUnexplainedEdges(unsigned long id) {
        unsigned ret = 0;
        const vector<unsigned> cluster = currentClusters[id].getElementsVector();
        for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
            vector<unsigned>::const_iterator it2 = it1;
            ++it2;
            for ( ; it2 != cluster.end(); ++it2) {
                if (!isThisEdgeExplained(*it1,*it2)) {
                    if (edgesToClusters[*it1][*it2] == 1) {
                        ++ret;
                    }
                }
            }
        }
        currentClusters[id].setUniquelyExplainedEdges(ret);
        return ret;
    }

    inline void setAllNumUniquelyExplained() {
        for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin();
             clustIt != currentClusters.end(); ++clustIt) {
            if ((clustIt->numElements() != 0) && (!clustIt->isValid())) {
                setNumUniquelyUnexplainedEdges(clustIt->getID());
            }
        }
    }


    inline void inactivateCluster(const unsigned long & clusterToInactivate, vector<string> & nodeIDsToNames, bool printClusterInfo = false) {
        /*quite similar to delete, but only need to remove from explanation, and set an inactive flag to that cluster*/
        if (currentClusters[clusterToInactivate].isAddedToExplain()) {
            removeClusterFromExplanation(clusterToInactivate);
        }
        currentClusters[clusterToInactivate].setInactive();
        if (printClusterInfo) {
            cout << "Inactivate: ";
           printCluster(currentClusters[clusterToInactivate].getElements(), nodeIDsToNames);
           cout << endl;
        }
    }

    inline void activateCluster(const unsigned long & clusterToActivate, vector<string> & nodeIDsToNames, bool printClusterInfo = false) {
        if (!currentClusters[clusterToActivate].isAddedToExplain()) {
            addClusterToExplanation(currentClusters[clusterToActivate].getElementsVector(), clusterToActivate);
        }
        currentClusters[clusterToActivate].setActive();
        if (printClusterInfo) {
            cout << "Activate: ";
            printCluster(currentClusters[clusterToActivate].getElements(), nodeIDsToNames);
            cout << endl;
        }
    }

    inline vector<unsigned long> deleteCluster(const unsigned long & clusterToDelete, vector<string> & nodeIDsToNames, bool printClusterInfo = true) {
        //TODO: should also delete edge from clusterGraph if this is the only cluster explaining that edge

        /*re-activate hidden cliques*/
        if (printClusterInfo) {
            cout << "Delete: ";
            printCluster(currentClusters[clusterToDelete].getElements(), nodeIDsToNames);
            cout << endl;
        }

        vector<unsigned long> hiddenClusters = currentClusters[clusterToDelete].getMergedFromID();
        vector<unsigned long> upperClusters = currentClusters[clusterToDelete].getMergedToID();

        if (currentClusters[clusterToDelete].isAddedToExplain()) {
            removeClusterFromExplanation(clusterToDelete);
        }
        openIDs.push_back(clusterToDelete);
        vector<ClusterBitset>::iterator clusterToDelete_it = currentClusters.begin();
        clusterToDelete_it += clusterToDelete;
        // it iterates through the nodes in the cluster to be deleted.  This allows us
        // to remove the now deleted cluster from the nodeToClusters structure
        for (unsigned i = 0; i < numNodes; ++i) {
            if (clusterToDelete_it->isElement(i)) {
                Utils::eraseInOrder(nodesToClusters[i], clusterToDelete);
            }
        }
        ++clustersDeleted;
        if (clusterToDelete_it->isNew()) {
            --newClusts;
        }
        *clusterToDelete_it = ClusterBitset();//this is needed

        /* remove the current cluster from the hidden clusters */
        for (vector<unsigned long>::iterator hiddenIt = hiddenClusters.begin(); hiddenIt !=  hiddenClusters.end(); ++hiddenIt) {
            currentClusters[*hiddenIt].removeMergedToID(clusterToDelete);
        }

        for (vector<unsigned long>::iterator upperIt = upperClusters.begin(); upperIt !=  upperClusters.end(); ++upperIt) {
            currentClusters[*upperIt].removeMergedFromID(clusterToDelete);
        }

        return hiddenClusters;
    }

    /*inline*/ void deleteClusters(vector<unsigned long> & clustersToDelete, vector<string> & nodeIDsToNames, graph_undirected_bitset & clusterGraph) {

        for (vector<unsigned long>::iterator clustersToDelete_it = clustersToDelete.begin();
             clustersToDelete_it != clustersToDelete.end(); ++clustersToDelete_it) {
            deleteCluster(*clustersToDelete_it, nodeIDsToNames);
        }

        return;
    }

    /*inline*/ void setClusterValid(const boost::dynamic_bitset<unsigned long> & cluster, graph_undirected_bitset & clusterGraph) {
        vector<unsigned> clusterElems;
        clusterElems.reserve(cluster.size());
        for (unsigned i = cluster.find_first(); i < cluster.size(); i = cluster.find_next(i)) {
            clusterElems.push_back(i);
        }
        for (unsigned i = 0; i < clusterElems.size()-1;  ++i) {
            for (unsigned j = i+1; j < clusterElems.size(); ++j) {
                setEdgeExplained(clusterElems[i], clusterElems[j]);
//                if (!clusterGraph.isEdge(clusterElems[i], clusterElems[j])) {//add edges to clusterGraph
//                   clusterGraph.addEdge(clusterElems[i], clusterElems[j]);
//                }
            }
        }
        return;
    }

    inline void setClusterValid(unsigned long clusterID, graph_undirected_bitset & clusterGraph) {
        setNumUniquelyUnexplainedEdges(clusterID);
        setClusterValid(getElements(clusterID), clusterGraph);
        currentClusters[clusterID].setValid();

        // don't deal with hidden clusters here. hidden clusters should be deleted alltogether in once
        return;
    }

    inline unsigned numNew() {
        return newClusts;
    }

    inline unsigned numClustersWithNode(unsigned nodeID) {
        return nodesToClusters[nodeID].size();
    }

    inline unsigned numCurrentClusters() {
        return currentClusters.size() - openIDs.size();
    }

    inline vector<unsigned long>::iterator clustersWithNodeBegin(unsigned nodeID) {
        return nodesToClusters[nodeID].begin();
    }

    inline vector<unsigned long>::iterator clustersWithNodeEnd(unsigned nodeID) {
        return nodesToClusters[nodeID].end();
    }


    /*inline*/ void sortNewClusters(vector<unsigned long> & sortedNewClusters) {
        /*sorting is still needed, to keep the familiar */
        vector<pair<unsigned long, unsigned> > newClustersAndCounts;
        newClustersAndCounts.reserve(numCurrentClusters());
        unsigned numFound = 0;
        for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
             clustIt != currentClusters.rend(); ++clustIt) {
            if ((clustIt->numElements() != 0) && (clustIt->isActive()) && (clustIt->isNew())) {
                newClustersAndCounts.push_back(make_pair(clustIt->getID(), clustIt->numElements()));
                ++numFound;
            }
        }//is this only new clusters?
        sort(newClustersAndCounts.begin(), newClustersAndCounts.end(), compPairSecondDescending);

        for (vector<pair<unsigned long, unsigned> >::iterator it = newClustersAndCounts.begin();
             it != newClustersAndCounts.end(); ++it) {
            sortedNewClusters.push_back(it->first);
        }
        return;
    }

    inline void clearEdgesToClusters() {
        for (unsigned i = 0; i < numNodes; ++i) {
            for (unsigned j = (i+1); j < numNodes; ++j) {
                edgesToClusters[i][j] = 0;
            }
        }
        return;
    }

    inline void addClustersToExplanations(vector<unsigned long> & sortedNewClusters) {
        for (vector<unsigned long>::iterator newClustIt = sortedNewClusters.begin();
             newClustIt != sortedNewClusters.end(); ++newClustIt) {
            if (currentClusters[*newClustIt].isActive() && !currentClusters[*newClustIt].isAddedToExplain()) {
                addClusterToExplanation(currentClusters[*newClustIt].getElementsVector(), *newClustIt);
            }
        }
        return;
    }


    /*inline*/ void prepareForValidityCheck(vector<unsigned long> & sortedNewClusters) {
        sortNewClusters(sortedNewClusters);
        addClustersToExplanations(sortedNewClusters);
        return;
    }

    bool isLargestExplainer(unsigned edgeEnd1, unsigned edgeEnd2, unsigned long id, vector<char> & idsChecked) {
        unsigned smallest = edgeEnd1;
        unsigned other = edgeEnd2;
        if (numClustersWithNode(edgeEnd2) < numClustersWithNode(edgeEnd1)) {
            smallest = edgeEnd2;
            other = edgeEnd1;
        }
        for (vector<unsigned long>::iterator nodeToClustersIt = clustersWithNodeBegin(smallest);
             nodeToClustersIt != clustersWithNodeEnd(smallest); ++nodeToClustersIt) {
            if (!idsChecked[*nodeToClustersIt] && currentClusters[*nodeToClustersIt].isElement(other)) {
                return false;
            }
        }
        return true;
    }

    inline unsigned getNumUnexplainedEdges(unsigned long id) {
        unsigned ret = 0;
        const vector<unsigned> cluster = currentClusters[id].getElementsVector();
        for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
            vector<unsigned>::const_iterator it2 = it1;
            ++it2;
            for ( ; it2 != cluster.end(); ++it2) {
                if (!isThisEdgeExplained(*it1,*it2)) {
                    //if (edgesToClusters[*it1][*it2] == 1) {
                    ++ret;
                    //}
                }
            }
        }
        return ret;
    }

    inline unsigned getNumUniquelyUnexplainedEdges(unsigned long id) {
        return currentClusters[id].getUniquelyExplainedEdges();
    }

    inline bool checkClusterValidity(unsigned long id) {
        const vector<unsigned> cluster = currentClusters[id].getElementsVector();
        //cout << id << ": ";
        for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
            vector<unsigned>::const_iterator it2 = it1;
            ++it2;
            for ( ; it2 != cluster.end(); ++it2) {
                //cout << edgesToClusters[*it1][*it2] << " " ;
                if (edgesToClusters[*it1][*it2] == 1) {// this edge only appears in one cluster. may use realedges so fewer necessary clusters?
//                    setNecessary(id);
                    return true;
                }
            }
        }
        return false;
    }

    inline const boost::dynamic_bitset<unsigned long> & getElements(unsigned long id) {
        return currentClusters[id].getElements();
    }

    inline unsigned numElements(unsigned long id) {
        return currentClusters[id].numElements();
    }

    inline void setOld(unsigned long id) {
        if (currentClusters[id].isNew()){
            --newClusts;
        }
        return currentClusters[id].setOld();
    }

    inline void setNew(unsigned long id) {
        if (!currentClusters[id].isNew()) {
            ++newClusts;
        }
        return currentClusters[id].setNew();
    }

    inline bool isNew(unsigned long id) {
        return currentClusters[id].isNew();
    }

    inline bool isActive(unsigned long id) {
        return currentClusters[id].isActive();
    }

    inline bool isValid(unsigned long id) {
        return currentClusters[id].isValid();
    }

    inline double getClusterWeight(unsigned long id) {
        return currentClusters[id].getWeight();
    }

    inline unsigned long maxClusterID() {
        return currentClusters.size();
    }

    inline double getCurWeight() {
        return curWeight;
    }

    inline double getMinWeightAdded() {
        return minWeightAdded;
    }

    inline double getThresh(unsigned long id) {
//        if (!currentClusters[id].isValid() && !currentClusters[id].isUnexplainedCounted()) {
//            setNumUniquelyUnexplainedEdges(id);
//        }
        return currentClusters[id].getThresh(alpha, getNumUnexplainedEdges(id));
    }

    inline double getMaxThresh() {
        unsigned numFound = 0;
        nextThreshold = 0;
        for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
             clustIt != currentClusters.rend(); ++clustIt) {
            if (clustIt->isNew()) {
                ++numFound;
                double clustThresh = getThresh(clustIt->getID());
                if (clustThresh > nextThreshold) {
                    nextThreshold = clustThresh;
                }
            }
            if (numFound == numNew()) {
                return nextThreshold;
            }
        }
        return nextThreshold;
    }

    void resetClusterWeight(unsigned long id, nodeDistanceObject & nodeDistances) {
        double weight = calculateClusterWeight(currentClusters[id].getElements(), nodeDistances);
        currentClusters[id].setWeight(weight);
        return;
    }

    unsigned numEdgesCovered() {
        graph_undirected_bitset coveredEdges(numNodes);
        for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
            if (clustIt->numElements() != 0) {
                vector<unsigned> clusterElems = clustIt->getElementsVector();
                for (unsigned i = 0; i < clusterElems.size()-1; ++i) {
                    for (unsigned j = i+1; j < clusterElems.size(); ++j) {
                        if (!coveredEdges.isEdge(clusterElems[i],clusterElems[j])) {
                            coveredEdges.addEdge(clusterElems[i],clusterElems[j]);
                        }
                    }
                }
            }
        }
        return coveredEdges.numEdges();
    }

    //TODO: this function is not good in some cases. Consider hide it
    bool isTooSmallForCurWeight(unsigned long id, unsigned lastLargestCluster, double lastLargestClusterWeight) {
        unsigned size = currentClusters[id].getElements().count();
        double latesmallThres_abs = ( log(numNodes) - log(size) ) /  log(numNodes) ;
        double latesmallThres_rel = lastLargestClusterWeight * ( log(numNodes) - log(size) ) / ( log(numNodes) - log(lastLargestCluster) ) ;
        double latesmallThres = min(latesmallThres_abs, latesmallThres_rel);
        if (curWeight < latesmallThres - 0.1) {
            return true;
        }
        else {
            return false;
        }
    }

    unsigned long getMaxClusterSizeForNode(unsigned long nodeID) {
        return maxClusterSizePerNodes[nodeID];
    }

    void updateMaxClusterSizeForNode(unsigned long nodeID, unsigned long size) {
        maxClusterSizePerNodes[nodeID] = size;
    }

    bool clusterSizeTooSmall(unsigned long id) {
        boost::dynamic_bitset<unsigned long> clust = currentClusters[id].getElements();
        unsigned long clustSize = clust.count();
        if (clustSize >= double(largestCluster) / 2) {
            return false; //don't need to look one by one
        }
        for (unsigned long i = clust.find_first(); i < clust.size(); i = clust.find_next(i)) {
            if (clustSize < double(getMaxClusterSizeForNode(i)) - 2) {
                return true;
            }
        }
        return false;
    }

    void setLargestClusterForNodes(unsigned long id) {
        boost::dynamic_bitset<unsigned long> clust = currentClusters[id].getElements();
        unsigned long clustSize = clust.count();
        for (unsigned long i = clust.find_first(); i < clust.size(); i = clust.find_next(i)) {
            if (clustSize > getMaxClusterSizeForNode(i)) {
                updateMaxClusterSizeForNode(i, clustSize);
            }
        }
    }

    vector<unsigned long> getMergeFromID(unsigned long id) {
        return currentClusters[id].getMergedFromID();
    }

    vector<unsigned long> getMergeToID(unsigned long id) {
        return currentClusters[id].getMergedToID();
    }

    void setMergedFromID(unsigned long id, vector <unsigned long> & mergeID) {
        currentClusters[id].setMergedFromID(mergeID);
    }

    void setMergedToID(unsigned long id, vector <unsigned long> & mergeID) {
        currentClusters[id].setMergeToID(mergeID);
    }

    void addMergedFromID(unsigned long id, unsigned long mergeID) {
        currentClusters[id].addMergedFromID(mergeID);
    }

    void addMergedToID(unsigned long id, unsigned long mergeID) {
        currentClusters[id].addMergedToID(mergeID);
    }

    void removeMergedFromID(unsigned long id, unsigned long mergeID) {
        currentClusters[id].removeMergedFromID(mergeID);
    }

    void removeMergedToID(unsigned long id, unsigned long mergeID) {
        currentClusters[id].removeMergedToID(mergeID);
    }

    unsigned clustersAdded;
    unsigned clustersDeleted;

private:
    vector<ClusterBitset> currentClusters;
    vector<unsigned long> openIDs;

    // Each node has an ID.  The nodesToClusters object can be accessed by that ID,
    // giving a set of iterators to cluster objects in the currentClusters list
    // which can then be derefenced to get a vector<unsigned> for the cluster
    vector<vector<unsigned long> > nodesToClusters;
    vector<vector<unsigned> > edgesToClusters;
    vector<vector<char> > isEdgeExplained;

    unsigned long numClusters;
    unsigned long nextID;
    unsigned newClusts;
    unsigned numNodes;
    double curWeight;
    double minWeightAdded;
    double firstWeight;
    double alpha;
    unsigned largestCluster;

    // Maximum weight of any new cluster curently in this list
    double maxNewWeight;
    double nextThreshold;

    // new feature: maximum cluster size per gene
    vector<unsigned long> maxClusterSizePerNodes;
};


namespace dagConstruct {

    /**/
    bool isClusterAncestor(const boost::dynamic_bitset<unsigned long> &cluster,
                           const boost::dynamic_bitset<unsigned long> &possibleAncestorCluster,
                           boost::dynamic_bitset<unsigned long> &unaccountedFor) {
        // Keep track of which genes in the ancestor cluster are "accounted for" (i.e. contained in) decsendent
        if (cluster.is_subset_of(possibleAncestorCluster)) {
            unaccountedFor -= cluster;
            return true;
        }
        return false;
    }

    bool isClusterAncestor(const vector<unsigned> &cluster, const vector<unsigned> &possibleAncestorCluster) {
        vector<unsigned>::const_iterator clusterIt = cluster.begin();
        vector<unsigned>::const_iterator ancestorIt = possibleAncestorCluster.begin();
        while ((clusterIt != cluster.end()) && (ancestorIt != possibleAncestorCluster.end())) {
            if (*ancestorIt < *clusterIt) {
                ++ancestorIt;
            } else if (*ancestorIt > *clusterIt) {
                return false;
            } else if (*ancestorIt == *clusterIt) {
                ++ancestorIt;
                ++clusterIt;
            }
        }
        if (clusterIt == cluster.end()) {
            return true;
        }
        return false;
    }



    inline unsigned performValidityCheck(currentClusterClassBitset & currentClusters,
                                     nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames) {
        vector<unsigned long> newClustersSorted; // where is this given?
        currentClusters.prepareForValidityCheck(newClustersSorted);
        unsigned deleted = 0;

        for (vector<unsigned long>::iterator newClusterIt = newClustersSorted.begin();
             newClusterIt != newClustersSorted.end(); ++newClusterIt) {
            if (!currentClusters.checkClusterValidity(*newClusterIt) && (currentClusters.isActive(*newClusterIt))) {
                currentClusters.deleteCluster(*newClusterIt, nodeIDsToNames, debug); // don't worry about the non-maximal child clusters. Since they won't pass here
                ++deleted;
            }
        }
        currentClusters.resetAllUnexplained(); // this only reset per cluster features

        return deleted;
    }

    bool isMinNodeDegreeMet(unsigned cluster1, unsigned cluster2, currentClusterClassBitset &currentClusters,
                            graph_undirected_bitset &clusterGraph, double density, vector<string> & nodeIDsToNames) {

        boost::dynamic_bitset<unsigned long> elements1 = currentClusters.getElements(cluster1);
        boost::dynamic_bitset<unsigned long> elements2 = currentClusters.getElements(cluster2);
        boost::dynamic_bitset<unsigned long> combination = currentClusters.getElements(cluster1);
        boost::dynamic_bitset<unsigned long> joint = currentClusters.getElements(cluster1);

        combination |= elements2;
        joint &= elements2;
        unsigned numCombined = combination.count();
        double numJoint = joint.count();
        double denom = numCombined - 1;
        double denom2 = numJoint - 1;

        // boundary condition can make this faster

        unsigned n1 = elements1.count();
        unsigned n2 = elements2.count();
        if (numJoint < 2) { return false;}

        //combine 1
//        for (unsigned i = elements1.find_first(); i < elements1.size(); i = elements1.find_next(i)) {
//            boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//            interactorsInCombo &= clusterGraph.getInteractors(i);
//            unsigned numInteractorsInCombo = interactorsInCombo.count();
//            if ((numInteractorsInCombo / denom) <= density) {
//                if (debug) { cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl; }
//                return false;
//            } else if (debug) { cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl; }
//        }
//        for (unsigned i = elements2.find_first(); i < elements2.size(); i = elements2.find_next(i)) {
//            if (elements1[i]) { continue; }
//            boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//            interactorsInCombo &= clusterGraph.getInteractors(i);
//            unsigned numInteractorsInCombo = interactorsInCombo.count();
//            if ((numInteractorsInCombo / denom) <= density) {
//                if (debug) { cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl; }
//                return false;
//            } else if (debug) { cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl; }
//        }

        //combine 2
        if ((numJoint > 1) && (numJoint/numCombined > density)) { // this faciliate merge of highly overlapped ones
            for (unsigned i = joint.find_first(); i < joint.size(); i = joint.find_next(i)) {
                boost::dynamic_bitset<unsigned long> interactorsInJoint = joint;
                interactorsInJoint &= clusterGraph.getInteractors(i);
                unsigned numInteractorsInJoint = interactorsInJoint.count();
                if ((numInteractorsInJoint / denom2) <= density ) {
                    if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInJoint / denom << endl;}
                    return false;
                }
                else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInJoint / denom << endl;}
            }
        }
        else if (n1 >= n2) {
            for (unsigned i = elements2.find_first(); i < elements2.size(); i = elements2.find_next(i) ) {
                if (joint[i]) { continue;}
                boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
                interactorsInCombo &= clusterGraph.getInteractors(i);
                unsigned numInteractorsInCombo = interactorsInCombo.count();
                if ((numInteractorsInCombo / denom) <= density ) {
                    if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
                    return false;
                }
                else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
            }
        }
        else {
            for (unsigned i = elements1.find_first(); i < elements1.size(); i = elements1.find_next(i) ) {
                if (joint[i]) { continue;}
                boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
                interactorsInCombo &= clusterGraph.getInteractors(i);
                unsigned numInteractorsInCombo = interactorsInCombo.count();
                if ((numInteractorsInCombo / denom) <= density ) {
                    if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
                    return false;
                }
                else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
            }
        }

        //combine 3


//        if (debug) {cout << n1 << "\t" << n2 << "\t" << denom << endl;}
//        if ( (n1-1)/denom > density) {
//            if (n1 >= n2) { // only test n2; otherwise, this means the above condition also applies for n2
//                for (unsigned i = elements2.find_first(); i < elements2.size(); i = elements2.find_next(i) ) {
//                    boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                    interactorsInCombo &= clusterGraph.getInteractors(i);
//                    unsigned numInteractorsInCombo = interactorsInCombo.count();
//                    if ((numInteractorsInCombo / denom) <= density ) {
//                        if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                        return false;
//                    }
//                    else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                }
//            }
//        }
//        else {
//            if ( (n2-1)/denom > density) {
//                //only test n1
//                for (unsigned i = elements1.find_first(); i < elements1.size(); i = elements1.find_next(i) ) {
//                    boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                    interactorsInCombo &= clusterGraph.getInteractors(i);
//                    unsigned numInteractorsInCombo = interactorsInCombo.count();
//                    if ((numInteractorsInCombo / denom) <= density ) {
//                        if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                        return false;
//                    }
//                    else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                }
//            }
//            else {//test both
//                if (n1<=n2) {//test the smaller one first
//                    for (unsigned i = elements1.find_first(); i < elements1.size(); i = elements1.find_next(i)) {
//                        boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                        interactorsInCombo &= clusterGraph.getInteractors(i);
//                        unsigned numInteractorsInCombo = interactorsInCombo.count();
//                        if ((numInteractorsInCombo / denom) <= density ) {
//                            if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                            return false;
//                        }
//                        else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                    }
//                    for (unsigned i = elements2.find_first(); i < elements2.size(); i = elements2.find_next(i)) {
//                        if (elements1[i]) { continue;}
//                        boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                        interactorsInCombo &= clusterGraph.getInteractors(i);
//                        unsigned numInteractorsInCombo = interactorsInCombo.count();
//                        if ((numInteractorsInCombo / denom) <= density ) {
//                            if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                            return false;
//                        }
//                        else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                    }
//                }
//                else {
//                    for (unsigned i = elements2.find_first(); i < elements2.size(); i = elements2.find_next(i)) {
//                        boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                        interactorsInCombo &= clusterGraph.getInteractors(i);
//                        unsigned numInteractorsInCombo = interactorsInCombo.count();
//                        if ((numInteractorsInCombo / denom) <= density ) {
//                            if (debug) {
//                                cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;
//                            }
//                            return false;
//                        }
//                        else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                    }
//                    for (unsigned i = elements1.find_first(); i < elements1.size(); i = elements1.find_next(i)) {
//                        if (elements2[i]) { continue;}
//                        boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
//                        interactorsInCombo &= clusterGraph.getInteractors(i);
//                        unsigned numInteractorsInCombo = interactorsInCombo.count();
//                        if ((numInteractorsInCombo / denom) <= density ) {
//                            if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                            return false;
//                        }
//                        else if (debug) {cout << nodeIDsToNames[i] << ' ' << numInteractorsInCombo / denom << endl;}
//                    }
//                }
//            }
//      }

        return true;
    }


    /*core clique finding algorithm. Don't know the detail*/
    void updateCliques(graph_undirected_bitset &clusterGraph,
                       vector<boost::dynamic_bitset<unsigned long> > &newClustersToAdd,
                       boost::dynamic_bitset<unsigned long> startClust, vector<unsigned> &neighborsOfBoth,
                       vector<char> &neighborsOfBothBits, unsigned i, vector<string> &nodeIDsToNames) {
        if (i == neighborsOfBoth.size()) {
            if (!Utils::insertInOrder(newClustersToAdd, startClust)) {
                cout << "# Repeated" << endl;
            }
            return;
        }
        if (startClust.is_subset_of(clusterGraph.getInteractors(neighborsOfBoth[i]))) {
            startClust[neighborsOfBoth[i]] = 1; //changed here. Why no ampersand?
            updateCliques(clusterGraph, newClustersToAdd, startClust, neighborsOfBoth, neighborsOfBothBits, i + 1,
                          nodeIDsToNames);
        } else {
            updateCliques(clusterGraph, newClustersToAdd, startClust, neighborsOfBoth, neighborsOfBothBits, i + 1,
                          nodeIDsToNames);

            vector<unsigned> T(clusterGraph.numNodes(), 0);

            boost::dynamic_bitset<unsigned long> nextClust = startClust;
            nextClust &= clusterGraph.getInteractors(neighborsOfBoth[i]);
            unsigned C_int_Ni_Count = nextClust.count();

            for (unsigned x = nextClust.find_first(); x < clusterGraph.numNodes(); x = nextClust.find_next(x)) {
                for (unsigned y = clusterGraph.getInteractors(x).find_first(); y < clusterGraph.numNodes();
                     y = clusterGraph.getInteractors(x).find_next(y)) {
                    if ((y != neighborsOfBoth[i]) && (startClust[y] == 0)) {
                        ++T[y];
                        if ((y < neighborsOfBoth[i]) && (T[y] == C_int_Ni_Count)) {
                            return;
                        }
                    }
                }
            }

            // Calculate S[y] = |N(y) intersection with (C - N(i))| for y not in C
            vector<unsigned> S(clusterGraph.numNodes(), 0);
            boost::dynamic_bitset<unsigned long> C_not_Ni = startClust;
            C_not_Ni -= clusterGraph.getInteractors(neighborsOfBoth[i]);
            unsigned C_not_Ni_count = C_not_Ni.count();
            for (unsigned x = C_not_Ni.find_first(); x < clusterGraph.numNodes(); x = C_not_Ni.find_next(x)) {
                for (unsigned y = clusterGraph.getInteractors(x).find_first(); y < clusterGraph.numNodes();
                     y = clusterGraph.getInteractors(x).find_next(y)) {
                    if ((startClust[y] == 0) && (neighborsOfBothBits[y] == 1)) {
                        ++S[y];
                    }
                }
            }

            unsigned k = 0;
            unsigned last_jk = 0;
            for (unsigned jk = C_not_Ni.find_first(); jk < clusterGraph.numNodes(); jk = C_not_Ni.find_next(jk)) {
                boost::dynamic_bitset<unsigned long> Njk_not_C = clusterGraph.getInteractors(jk);
                for (unsigned y = Njk_not_C.find_first(); y < neighborsOfBoth[i]; y = Njk_not_C.find_next(y)) {
                    if ((T[y] == C_int_Ni_Count) && (neighborsOfBothBits[y] == 1)) {
                        if (y >= jk) {
                            --S[y];

                        } else if (((S[y] + k) == C_not_Ni_count) && (y >= last_jk)) {
                            return;
                        }
                    }
                }
                last_jk = jk;
                ++k;
            }

            unsigned jp = last_jk;
            if (jp < (neighborsOfBoth[i] - 1)) {
                return;
            }
            unsigned inNext = nextClust.find_first();
            for (unsigned y = clusterGraph.getInteractors(inNext).find_first(); y < clusterGraph.numNodes();
                 y = clusterGraph.getInteractors(inNext).find_next(y)) {
                if ((y < neighborsOfBoth[i]) && (startClust[y] == 0) && (T[y] == C_int_Ni_Count) && (S[y] == 0) &&
                    (jp < y)) {
                    return;
                }
            }
            nextClust[neighborsOfBoth[i]] = 1;
            updateCliques(clusterGraph, newClustersToAdd, nextClust, neighborsOfBoth, neighborsOfBothBits, i + 1,
                          nodeIDsToNames);
        }
    }


    /**/
    bool findClustsWithSeed(boost::dynamic_bitset<unsigned long> seedClust,
                            graph_undirected_bitset &clusterGraph,
                            vector<boost::dynamic_bitset<unsigned long> > &newClustersToAdd,
                            vector<string> &nodeIDsToNames) {

        // fast filter of the edges that are already in a new cluster
        for (unsigned i = 0; i < newClustersToAdd.size(); ++i) {
            if (seedClust.is_subset_of(newClustersToAdd[i])) {
                return true;
            }
        }

        // to decide if a node is the neighbor of all nodes in seedClust. If true, then this node needs to be merge into seedClust
        vector<unsigned> neighborsOfAll;
        vector<char> neighborsOfAllBits(clusterGraph.numNodes(), 0);
        neighborsOfAll.reserve(clusterGraph.numNodes());


        for (unsigned i = 0; i < clusterGraph.numNodes(); ++i) {// I think here I can only consider affected nodes
            if ((seedClust[i] == 0) && (seedClust.is_subset_of(clusterGraph.getInteractors(i)))) {
                neighborsOfAll.push_back(i);
                neighborsOfAllBits[i] = 1;
            }
        }

        if (neighborsOfAll.size() > 0) {
            updateCliques(clusterGraph, newClustersToAdd, seedClust, neighborsOfAll, neighborsOfAllBits, 0,
                          nodeIDsToNames); //what is changed?
            return true;
        }

        return false;
    }


    /**/
    void updateClustersWithEdges(vector<pair<unsigned, unsigned> > &edgesToAdd,
                                     currentClusterClassBitset &currentClusters, graph_undirected_bitset &clusterGraph,
                                     nodeDistanceObject &nodeDistances, vector<string> &nodeIDsToNames, bool withRealEdges = true) {

        cout << "# Adding " << edgesToAdd.size() << " edges" << endl;

        // add edges, and check affected nodes
        vector<char> affectedNodes(clusterGraph.numNodes(), 0);
        for (unsigned long edgesToAddCounter = 0; edgesToAddCounter != edgesToAdd.size(); ++edgesToAddCounter) {

            if (!clusterGraph.isEdge(edgesToAdd[edgesToAddCounter].first,
                                     edgesToAdd[edgesToAddCounter].second)) {
                clusterGraph.addEdge(edgesToAdd[edgesToAddCounter].first,
                                     edgesToAdd[edgesToAddCounter].second);
            }
            affectedNodes[edgesToAdd[edgesToAddCounter].first] = 1;
            affectedNodes[edgesToAdd[edgesToAddCounter].second] = 1;
        }

        unsigned long numAff = 0;
        for (vector<char>::iterator thisIt = affectedNodes.begin(); thisIt != affectedNodes.end(); ++thisIt) {
            if (*thisIt) {
                ++numAff;
            }
        }
        cout << "# Nodes affected = " << numAff << endl;

        vector<boost::dynamic_bitset<unsigned long> > newClustersToAdd;
        newClustersToAdd.reserve(20000);
        vector<char> affectedClusters(currentClusters.maxClusterID(), 0);

        // default behavior (no chordal graph)

        // label affected clusters
        for (unsigned i = 0; i < affectedNodes.size(); ++i) {
            if (affectedNodes[i]) {
                for (vector<unsigned long>::iterator cNodesIt = currentClusters.clustersWithNodeBegin(i);
                     cNodesIt != currentClusters.clustersWithNodeEnd(i); ++cNodesIt) {
                    affectedClusters[*cNodesIt] = 1;
                }
            }
        }

        numAff = 0;
        for (vector<char>::iterator thisIt = affectedClusters.begin();
             thisIt != affectedClusters.end(); ++thisIt) {
            if (*thisIt) {
                ++numAff;
            }
        }
        cout << "# Clusters affected = " << numAff << endl;

        unsigned inactivated = 0;

        // expanding existing clusters, in currentClusters
        for (unsigned i = 0; i < affectedClusters.size(); ++i) {
            if (affectedClusters[i] && currentClusters.isNew(i) && currentClusters.isActive(i)) {
                if (findClustsWithSeed(currentClusters.getElements(i), clusterGraph, newClustersToAdd,
                                       nodeIDsToNames)) {
                    //TODO: change this to cluster inactivation
                    currentClusters.inactivateCluster(i, nodeIDsToNames, debug); //set true for now for debugging
                    ++inactivated;
                }
            }
        }
        cout << "# Num of non-maximal clusters inactivated here: " << inactivated << endl;

        // adding brand new clusters, using new edges added this round
        for (unsigned long edgesToAddCounter = 0; edgesToAddCounter != edgesToAdd.size(); ++edgesToAddCounter) {
            boost::dynamic_bitset<unsigned long> seedClust(clusterGraph.numNodes());
            seedClust[edgesToAdd[edgesToAddCounter].first] = 1;
            seedClust[edgesToAdd[edgesToAddCounter].second] = 1;
            if (!findClustsWithSeed(seedClust, clusterGraph, newClustersToAdd,
                                    nodeIDsToNames)) { //seedClust is updated inside the findClustWithSeed function
                // seedClust is a maximal clique of just two nodes
                Utils::insertInOrder(newClustersToAdd, seedClust);
            }
        }

        /* add clusters */
        cout << "# Found " << newClustersToAdd.size() << " new clusters to add" << endl;
        for (vector<boost::dynamic_bitset<unsigned long> >::iterator clustersToAddIt = newClustersToAdd.begin();
             clustersToAddIt != newClustersToAdd.end(); ++clustersToAddIt) {
//                cout << "# adding Clusters" << endl;
//                printCluster(*clustersToAddIt, nodeIDsToNames);
//                cout << endl;
            currentClusters.addCluster(*clustersToAddIt, nodeDistances, nodeIDsToNames);
        }
        //done
    }

    bool combineClusters(
            vector<pair<pair<unsigned long, unsigned long>, double>> &clustersToCombine,
            currentClusterClassBitset &currentClusters,
            graph_undirected_bitset & clusterGraph,
            vector<string> &nodeIDsToNames,
            nodeDistanceObject &nodeDistances) {

        vector<pair<unsigned, unsigned> > edgesToAdd;
        edgesToAdd.reserve(clusterGraph.numEdges());

        vector<vector<char> > edgeWillBeAdded(nodeDistances.numNodes(), vector<char>(nodeDistances.numNodes(), false));
        sort(clustersToCombine.begin(), clustersToCombine.end(), compClustersToCombine); // not sure whether this sort is helpful

        for (vector<pair<pair<unsigned long, unsigned long>, double>>::iterator clustersToCombineIt = clustersToCombine.begin(); // why is this complaining
             clustersToCombineIt != clustersToCombine.end(); ++clustersToCombineIt) {

            boost::dynamic_bitset<unsigned long> onlyIn1 = currentClusters.getElements(
                    clustersToCombineIt->first.first);
            boost::dynamic_bitset<unsigned long> onlyIn2 = currentClusters.getElements(
                    clustersToCombineIt->first.second);
            onlyIn1 -= currentClusters.getElements(clustersToCombineIt->first.second);
            onlyIn2 -= currentClusters.getElements(clustersToCombineIt->first.first);

            for (unsigned i = onlyIn1.find_first(); i < onlyIn1.size(); i = onlyIn1.find_next(i)) {
                boost::dynamic_bitset<unsigned long> newEdgesFrom1To2 = onlyIn2;
                newEdgesFrom1To2 -= clusterGraph.getInteractors(
                        i); //excluding old edges (THINK ABOUT IT IS clusterGraph or clusterGraph)
                for (unsigned j = newEdgesFrom1To2.find_first();
                     j < newEdgesFrom1To2.size(); j = newEdgesFrom1To2.find_next(j)) {
                    if (!edgeWillBeAdded[i][j]) { // don't comment out this line
                        edgesToAdd.push_back(make_pair(i, j)); //here you are adding all edges!
                        edgeWillBeAdded[i][j] = true;
                        edgeWillBeAdded[j][i] = true;
                    }
                }
            }
        }

        if (edgesToAdd.size() == 0) {
            return false;
        }

        updateClustersWithEdges(edgesToAdd, currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);

        return true;
    }

    /**/
    bool addMissingEdges(currentClusterClassBitset &currentClusters,
                         double density, double threshold,
                         graph_undirected_bitset & clusterGraph,
                         vector<string> &nodeIDsToNames,
                         nodeDistanceObject &nodeDistances,
                         graph_undirected_bitset &realEdges) {

        vector<pair<pair<unsigned long, unsigned long>, double>> clustersToCombine; //the bitset is the element after combine

        unsigned long maxClusterID = currentClusters.maxClusterID();
        vector<bool> clustersChecked(maxClusterID, false);

//        int nvalidCluster = 0;
//        cout << "debug here" << endl;
//        cout << currentClusters.numCurrentClusters() << endl;
//        cout << maxClusterID << endl;
        for (unsigned long i = 0; i < maxClusterID; ++i) {
            if (currentClusters.isNew(i) &&
                    currentClusters.isActive(i) &&
                    (currentClusters.numElements(i) !=0) ){//&& //don't exactly why I need this,  seems to be duplicated to 1,2
//                    (currentClusters.getThresh(i) >= currentClusters.getCurWeight()))  {
                clustersChecked[i] = true;
                for (unsigned long j = 0; j < maxClusterID; ++j) {
                    /*I think I should consider all clusters*/
                    if ((j != i) && (!clustersChecked[j]) &&
                        currentClusters.isActive(j) &&
                        (currentClusters.numElements(j) != 0)) {//&&
//                            (currentClusters.getThresh(j) >= currentClusters.getCurWeight() ))  {
//                        if (currentClusters.isValid(i) || currentClusters.isValid(j)) {
//                            nvalidCluster = nvalidCluster + 1;
//                        }
//                        boost::dynamic_bitset<unsigned long> proposedCombinedCluster(realEdges.numNodes(), 0);
                        if (debug) {
                            cout << "# Merging attempt" << i << " " << j << endl;
                            printCluster(currentClusters.getElements(i), nodeIDsToNames);
                            cout << "\t";
                            printCluster(currentClusters.getElements(j), nodeIDsToNames);
                            cout << endl;

                        }
                        if (isMinNodeDegreeMet(i, j, currentClusters, realEdges, density,
                                               nodeIDsToNames)) {// try clusterGraph in this new version
                            double tempWeight;
                            if (debug) {
//                                printCluster(currentClusters.getElements(i), nodeIDsToNames);
//                                cout << "\t";
//                                printCluster(currentClusters.getElements(j), nodeIDsToNames);
//                                cout << endl;
                                cout << "# Successful merging" << endl;
                            }
                            if (currentClusters.getThresh(j) > currentClusters.getThresh(i)) {
                                tempWeight = currentClusters.getThresh(i);
                            } else {
                                tempWeight = currentClusters.getThresh(j);
                            }
                            clustersToCombine.push_back(make_pair(make_pair(i, j), tempWeight));
                        }
                    }
                }
            }
        }
//        if (debug && (nvalidCluster > 1)) {cout << "# valid cluster is considered" << endl;}
        cout << "# Considering combining " << clustersToCombine.size() << " pairs of clusters" << endl;
        if (clustersToCombine.size() > 0) {
//            cout << "# Current Weight is " << curWeight << endl;
            bool realCombine = combineClusters(clustersToCombine, currentClusters, clusterGraph, nodeIDsToNames, nodeDistances);
            return realCombine; // if return true, new cluster are created. Should consider more combine
        }
        return false;
    }


    // ** Main function ** //
    unsigned getFuzzyClustersBitsetThreshold(nodeDistanceObject &nodeDistances, map<string, unsigned> &nodeNamesToIDs,
                                             vector<validClusterBitset> & validClusters, double alpha = 0,
                                             double beta = 1) {

        // * profiling * //
        time_t start, end;
        time(&start);
        double dif;

        // track
        unsigned numRealEdgesAdded = 0;
        unsigned numRealEdgesLastRound = 0;
        unsigned largestCluster = 0;
        unsigned lastLargestCluster = 0;

        // * node (gene) names * //
        vector<string> nodeIDsToNames(nodeNamesToIDs.size(), string(""));
        for (map<string, unsigned>::iterator it = nodeNamesToIDs.begin(); it != nodeNamesToIDs.end(); ++it) {// what is going on here
            nodeIDsToNames[it->second] = it->first;
        }
        unsigned numNodes = nodeDistances.numNodes();

        // * define graph objects * //clusterGraphlastRoundEdges
        graph_undirected_bitset clusterGraph(numNodes); // all the edges covered by clusters by this time
        graph_undirected_bitset realEdges(numNodes); // all the real edges by this time

        //
        double dt = nodeDistances.sortedDistancesBegin()->second; // Current threshold, starting with the maximum similarity
        double last_dt = dt;
        double lastLargestClusterWeight = dt;
        currentClusterClassBitset currentClusters(numNodes, dt, alpha);

        sortedDistanceStruct::iterator distanceIt = nodeDistances.sortedDistancesBegin();
        unsigned totalEdges = numNodes * (numNodes - 1) / 2;

        // iterative clique finding
        while ((realEdges.numEdges() < 0.5 * totalEdges) &&
               (distanceIt != nodeDistances.sortedDistancesEnd())) { // termination conditions

            unsigned numRealEdgesThisRound = 0;
//            clusterGraph = realEdges; // is this making a copy

            vector<pair<unsigned, unsigned> > edgesToAdd;
            edgesToAdd.reserve(20000); //don't know if this is helpful

            double addUntil = dt;
//            while (numRealEdgesThisRound <= numRealEdgesLastRound) {
            addUntil = dt - alpha; // addUntil can cross multiple alpha if there is a region with low edge density
            while ((distanceIt != nodeDistances.sortedDistancesEnd()) &&
                   (distanceIt->second >= addUntil)) {
                unsigned firstNode = distanceIt->first.first;
                unsigned secondNode = distanceIt->first.second;
                realEdges.addEdge(firstNode,
                                  secondNode); // add edges to realEdges; this is the only place that realEdges change
                ++numRealEdgesAdded;
                ++numRealEdgesThisRound;
                edgesToAdd.push_back(make_pair(firstNode, secondNode));
                ++distanceIt;
            }
            last_dt = dt;
            if (distanceIt != nodeDistances.sortedDistancesEnd()) {
                dt = addUntil; //dt is already moved to the next level
            } else {
                dt = 0;
                break;
            }
//            }


            cout << "# Current distance: " << distanceIt->second << "\t" << "Add until: " << addUntil << "\t" << endl;
            cout << "# Num of real edges added: " << numRealEdgesAdded << endl;
//            cout << "# Use real edges to expand clusters" << endl;
            updateClustersWithEdges(edgesToAdd, currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);

            time(&end);
            dif = difftime(end, start);
            cout << "# Time elapsed: " << dif << " seconds" << endl;

            numRealEdgesLastRound = numRealEdgesThisRound;

            currentClusters.setCurWeight(dt);


            unsigned long numClustersBeforeDelete = currentClusters.numCurrentClusters();
            cout << "# Number of clusters before merging: " << numClustersBeforeDelete << endl;

            vector<unsigned long> newClustersSorted;

            if (beta < 1) { // think about chordal later

                cout << "# Adding missing edges...checking " << currentClusters.numCurrentClusters() << " cliques" << endl;
                //TODO: should only consider active clusters here
                unsigned ndeleted = performValidityCheck(currentClusters, nodeDistances, nodeIDsToNames);
                bool newEdgesAdded = addMissingEdges(currentClusters, beta, alpha, clusterGraph, nodeIDsToNames, nodeDistances, realEdges);

                cout << "# " << ndeleted << " clusters deleted in validity check; Now has " << currentClusters.numCurrentClusters() << " clusters" << endl;
                time(&end);
                dif = difftime(end, start);
                cout << "# Time elapsed: " << dif << " seconds" << endl;

                while (newEdgesAdded == true) { //recursively merging cliques
                    // delete invalid clusters after merging. remember, if a merged clique is deleted, the two original cliques must be returned
                    unsigned ndeleted = performValidityCheck(currentClusters, nodeDistances, nodeIDsToNames);
                    newEdgesAdded = addMissingEdges(currentClusters, beta, alpha, clusterGraph, nodeIDsToNames, nodeDistances, realEdges);
                    cout << "# " << ndeleted << " clusters deleted in validity check; Now has  " << currentClusters.numCurrentClusters() << " clusters" << endl;
                }
            }
            currentClusters.sortNewClusters(newClustersSorted); //sort by size, ascending; this only has new active clusters
//                currentClusters.prepareForValidityCheck(newClustersSorted); // sort and add edge explaination

            cout << "# Current number of clusters:" << currentClusters.numCurrentClusters() << endl;
            cout << "# New clusters to evaluate:" << newClustersSorted.size() << endl;
            time(&end);
            dif = difftime(end, start);
            cout << "# Time elapsed: " << dif << " seconds" << endl;


            vector<unsigned long> tempNewAndValid;

//            unsigned weight_filter = 0;
            unsigned small_filter = 0;
            unsigned unique_filter = 0;

            unsigned n_pass_filter = 0;

            while (newClustersSorted.size()>0) {
//                while ((newClustersSorted.size() > 0) &&
//                        (currentClusters.getMaxThresh() >= dt)) { //TODO: this explains why the first a few iteration is wrong
//                    currentClusters.sortNewClusters(newClustersSorted);
                unsigned long clusterTop = newClustersSorted.back(); // think about the order
                newClustersSorted.pop_back();
                if (currentClusters.numElements(clusterTop) ==0 ) { //some clusters are "hold-on" this stage.
                    continue;
                }

                //TODO: this function can go to updateClusterWeight (downside: prohitbit merging. left it as is for now)
                /*filter : see if the term is too small for the current weight*/

//                if (currentClusters.isTooSmallForCurWeight(clusterTop, lastLargestCluster, lastLargestClusterWeight)) {
                if (currentClusters.clusterSizeTooSmall(clusterTop)) {
                    if (debug) {
                        cout << "Size too small: ";
                    }
                    currentClusters.deleteCluster(clusterTop, nodeIDsToNames, debug);

                    //don't have to delete hiddenClusters here due to the final step
                    ++small_filter;
                    continue;
                }

                // weight filter
//                if (currentClusters.getThresh(clusterTop) < dt) {
//                    //if weight is not qualified (this happen often when merge happens; also delete and break down
//                    if (debug) {
//                        cout << "Weight failed: ";
//                    }
//                    vector<unsigned long> hiddenClusters = currentClusters.deleteCluster(clusterTop, nodeIDsToNames, debug);
//                    for (vector<unsigned long>::iterator hidden_it = hiddenClusters.begin();
//                         hidden_it != hiddenClusters.end(); ++hidden_it) {
//                        if ((currentClusters.isNew(*hidden_it)) &&
//                                (currentClusters.getMergeToID(*hidden_it).size() == 0)){
//                            //make it active again
//                            if (debug) {
//                                cout << "Activate again: ";
//                            }
//                            currentClusters.activateCluster(*hidden_it, nodeIDsToNames, debug);
//                            newClustersSorted.push_back(*hidden_it);
//                        }
//                    }
//                    ++weight_filter;
//                    continue;
//                }

                /*see if the term does not have many unique edges*/

                currentClusters.setNumUniquelyUnexplainedEdges(clusterTop);
//                    double uniqueThresh = (1-beta) * currentClusters.getElements(clusterTop).count();
                double uniqueThresh =  0.5 * (currentClusters.getElements(clusterTop).count() - 1) + 0.05* pow(currentClusters.getElements(clusterTop).count(), 2.0); //soft for small, I think quite strong for big

                if (currentClusters.getNumUniquelyUnexplainedEdges(clusterTop) < uniqueThresh) {
                    if (debug) {
                        cout << "Uniqueness failed: ";
                        cout << currentClusters.getElements(clusterTop).count() << "\t" << currentClusters.getNumUniquelyUnexplainedEdges(clusterTop) << "\t" << uniqueThresh
                             << "\t" << endl;
                    }
                    vector<unsigned long> hiddenClusters = currentClusters.deleteCluster(clusterTop, nodeIDsToNames, debug);
//                    for (vector<unsigned long>::iterator hidden_it = hiddenClusters.begin();
//                         hidden_it != hiddenClusters.end(); ++hidden_it) {
//                        if ((currentClusters.isNew(*hidden_it)) &&
//                                (currentClusters.getMergeToID(*hidden_it).size() == 0)){
//                            //make it active again
//                            if (debug) {
//                                cout << "Activate again: ";
//                            }
//                            currentClusters.activateCluster(*hidden_it, nodeIDsToNames, debug);
//                            newClustersSorted.push_back(*hidden_it);
//                        }
//                    }
                    ++unique_filter;
                    continue;
                }
                /*note that if a term is deleted, need trace back and recheck some smaller terms*/
                /*pass all filter*/
                ++n_pass_filter;
//                    tempNewAndValid.push_back(clusterTop);
            }

//            cout << "# " << weight_filter << " clusters failed the weight filter" << endl;
            cout << "# " << small_filter << " clusters failed the late-and-small filter" << endl;
            cout << "# " << unique_filter << " clusters failed the uniqueness filter" << endl;
            cout << "# " << n_pass_filter << " clusters are valid" << endl;
            //sort again by size
            currentClusters.sortNewClusters(tempNewAndValid);

            unsigned numNonMaximalDeleted = 0;
            for (vector<unsigned long>::iterator newValidClusterIt = tempNewAndValid.begin();
                 newValidClusterIt != tempNewAndValid.end(); ++newValidClusterIt) {
                double clustWeight = currentClusters.getClusterWeight(*newValidClusterIt);

                //treat 2 and 3 genes terms differently
                if (currentClusters.getElements(*newValidClusterIt).count() < 4) {
                    if (clustWeight < dt + alpha) {
                        continue;
                    }
                }

                validClusters.push_back(
                        validClusterBitset(currentClusters.getElements(*newValidClusterIt), 0, clustWeight));
                cout << "# Valid cluster:\t";
                printCluster(currentClusters.getElements(*newValidClusterIt), nodeIDsToNames);
                currentClusters.setClusterValid(*newValidClusterIt, realEdges); //successful clusters should be seen as real edges in future (okay because their weights doesn't change)
                cout << "\t" << clustWeight << "\t"
                     << currentClusters.getNumUniquelyUnexplainedEdges(*newValidClusterIt) << "\t" << last_dt
                     << endl;
                // now it is safe to delete the hidden cluster of the valid cluster. see the change in setClusterValid.
                if (validClusters.back().numElements() > largestCluster) {
                    lastLargestCluster = largestCluster;
                    lastLargestClusterWeight = clustWeight;
                    largestCluster = validClusters.back().numElements();
                    currentClusters.setLargestCluster(largestCluster);
                }
                currentClusters.setOld(*newValidClusterIt);
                currentClusters.setLargestClusterForNodes(*newValidClusterIt);

                vector<unsigned long> nonMaximals = currentClusters.getMergeFromID(*newValidClusterIt);
                for (vector<unsigned long>::iterator nonMaximalIt = nonMaximals.begin();
                 nonMaximalIt != nonMaximals.end(); ++nonMaximalIt) {
                    if ( !currentClusters.isValid(*nonMaximalIt) ) {
                        currentClusters.deleteCluster(*nonMaximalIt, nodeIDsToNames, debug);
                    }
                    ++numNonMaximalDeleted;
                }

            }
            cout << "# Non-maximal finally deleted: " << numNonMaximalDeleted << endl;


            cout << "# dt: " << last_dt << endl;
            cout << "# Next dt: " << dt << endl;
//            cout << "# Num current clusters before delete: " << numClustersBeforeDelete << endl;
            cout << "# Num current clusters: " << currentClusters.numCurrentClusters() << endl;
            cout << "# Num valid clusters: " << validClusters.size() << endl;
            cout << "# Largest cluster: " << largestCluster << endl;
            cout << "# Num edges in clusterGraph: " << clusterGraph.numEdges() << endl;
            cout << "# Num real edges added: " << numRealEdgesAdded << endl;
            cout << "# Num edges covered by valid clusters: " << realEdges.numEdges() << endl;
//                cout << "# Num edges inferred: " << clusterGraph.numEdges() - clusterGraph.numEdges() << endl;
            time(&end);
            dif = difftime(end, start);
            cout << "# Time elapsed: " << dif << " seconds" << endl;
        }
        return currentClusters.numCurrentClusters();
    }


    // Take list of valid clusters and turn it into an ontology, assuming the perfect case
    // where all edges in a cluster are identical
    void clustersToDAG(vector<validClusterBitset> validClusters, DAGraph & ontology, unsigned numTerminalNodes) {
        //unsigned numTerminalNodes = nodeDistances.numNodes();

        sort(validClusters.begin(), validClusters.end());

        ontology.reserveNodes(ontology.numNodes() + validClusters.size() + 1);
        unsigned clustCount = 0;
        vector<validClusterBitset>::iterator clustersIt = validClusters.begin();
        unsigned firstNodeID = ontology.addNode();
        clustersIt->setID(firstNodeID);
        for (unsigned i = 0; i < numTerminalNodes; ++i) {
            if (clustersIt->isElement(i)) {
                ontology.addEdge(firstNodeID, i);
            }
        }
        ontology.setWeight(clustersIt->getWeight(),firstNodeID);
        double geneWeight = clustersIt->getWeight();
        ++clustersIt;
        ++clustCount;
        for ( ; clustersIt != validClusters.end(); ++clustersIt) {
            //cout << clustCount << endl;
            unsigned newNodeID = ontology.addNode();
            clustersIt->setID(newNodeID);

            boost::dynamic_bitset<unsigned long> unaccountedFor = clustersIt->getElements();

            for (vector<validClusterBitset>::reverse_iterator possibleDescendentsIt(clustersIt);
                 possibleDescendentsIt != validClusters.rend(); ++possibleDescendentsIt) {
                if (possibleDescendentsIt->numElements() < clustersIt->numElements()) {
                    // A cluster is a child if it is not already a
                    // descendent of the current cluster via some other cluster,
                    // and all of its genes are contained in the possible ancestor cluster,

                    bool isAlreadyDescendent = ontology.isDescendent(possibleDescendentsIt->getID(), newNodeID);
                    if ((!isAlreadyDescendent) &&
                        (isClusterAncestor(possibleDescendentsIt->getElements(), clustersIt->getElements(), unaccountedFor))) {
                        //cout << "adding edge " << newNodeID << "\t" << possibleDescendentsIt->getID() << endl;
                        ontology.addEdge(newNodeID, possibleDescendentsIt->getID());
                    }
                }
            }

            // Add edges to genes which are not contained in any of the child nodes
            for (unsigned i = 0; i < numTerminalNodes; ++i) {
                if (unaccountedFor[i] == 1) {
                    ontology.addEdge(newNodeID, i);
                }
            }

            // Set weight of cluster
            ontology.setWeight(clustersIt->getWeight(),newNodeID);
            if (geneWeight < clustersIt->getWeight()) {
                geneWeight = clustersIt->getWeight();
            }
            //cout << ontology.getName(newNodeID) << "\t" << ontology.getWeight(newNodeID) << endl;
            ++clustCount;
        }

        if (geneWeight > 1) {
            ontology.setGeneWeight(geneWeight);
        } else {
            ontology.setGeneWeight(1);
        }

        //cout << "Finished creating DAG - adding root node" << endl;
        unsigned firstTopLevelNode;
        bool firstTopLevelFound = false;
        bool secondTopLevelFound = false;
        long rootID = -1;

        unsigned numTopLevel = 0;
        for (vector<DAGNode>::iterator nodeIt = ontology.nodesBegin(); nodeIt != ontology.nodesEnd(); ++nodeIt) {
            if ((nodeIt->numParents() == 0) && (nodeIt->getID() != rootID)) {
                //cout << nodeIt->getID() << " " << nodeIt->getName() << " is top level" << endl;
                ++numTopLevel;
                //cout << numTopLevel << endl;
                if (!firstTopLevelFound) {
                    firstTopLevelFound = true;
                    firstTopLevelNode = nodeIt->getID();
                } else if (!secondTopLevelFound) {
                    secondTopLevelFound = true;
                    unsigned curItPos = nodeIt->getID();

                    rootID = ontology.addNode();
                    // Adding node may invalidate nodeIt, so fix it
                    nodeIt = ontology.nodesBegin();
                    nodeIt += curItPos;

                    ontology.setWeight(0, rootID);
                    ontology.addEdge(rootID, firstTopLevelNode);
                    ontology.addEdge(rootID, nodeIt->getID());
                } else {
                    if (rootID != nodeIt->getID()) {
                        ontology.addEdge(rootID, nodeIt->getID());
                    }
                }
            }
        }
        return;
    }

    // Main clustering function - cluster input graph into Ontology
    void constructDAG(graph_undirected & input_graph, DAGraph & ontology, nodeDistanceObject & nodeDistances, double threshold, double density, bool debugging) {
        //ontology = DAGraph();
        map<string,unsigned> geneNamesToIDs;

        if (debugging) {
            debug = true;
        }

        // Add all nodes in the input network to the ontology as genes
        for (vector<Node>::iterator nodesIt = input_graph.nodesBegin(); nodesIt != input_graph.nodesEnd(); ++nodesIt) {
            geneNamesToIDs[nodesIt->getName()] = nodesIt->getID();
            ontology.addNode(nodesIt->getName(), geneNamesToIDs);
        }

        nodeDistances = nodeDistanceObject(input_graph);

        vector<validClusterBitset> validClusters;
        dagConstruct::getFuzzyClustersBitsetThreshold(nodeDistances, geneNamesToIDs, validClusters, threshold, density); //all important stuff in here

        // Got the clusters - build ontology assuming any clusters wholly contained in others are descendents
        cout << "# Num clusters: " << validClusters.size() << endl;
        dagConstruct::clustersToDAG(validClusters, ontology, geneNamesToIDs.size());
        return;
    }
}


#endif //CLIXO_0_3_DAGCONSTRUCT_H
