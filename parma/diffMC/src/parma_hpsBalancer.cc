#include <apfPartition.h>
#include <PCU.h>
#include "parma_sides.h"
#include "parma_entWeights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "zeroOneKnapsack.h"
#include "maximalIndependentSet/mis.h"
#include <limits>

using std::vector;
using std::set;

namespace parma {
  //header for avgWeights
  double avgWeight(Weights* w);

  class MergeTargets {  // we don't really need a map/associative container here - a list/vector/array would work
    public:
      //Have this storing the results in Targets associative class and assuming maxW = avgWeight * maxImb
      //maxW is also = HeavyImb
      MergeTargets(Sides* s, Weights* w, double maxW)
      {
        //If part is heavy or empty, exit function
        if (w->self() >= maxW || w->self() == 0) return;
        PCU_Debug_Print("Part %d of weight %f is light and not empty with imb of %f, %u neighbors--\n", PCU_Comm_Self(), w->self(), maxW, w->size());
        
        int* nborPartIds = new int[w->size()];
        int i = 0;
        const Weights::Item* weight;        
        w->begin();         
        while( (weight = w->iterate()) ) 
          nborPartIds[i++] = weight->first;    
        w->end();         

        //iterating through the neighbor weights to determine the minimum weight
        double minWeight = std::numeric_limits<double>::max(); 
        w->begin(); 
        while( (weight = w->iterate()) ) 
          if ( weight->second < minWeight )
            minWeight = weight->second;
        w->end();
        //normalizing the neighbor weights with a dividing factor to increase the accuracy of knapsack
        double t0 = MPI_Wtime();
        bool divide = true;
        int* normalizedIntWeights = new int[w->size()];
        unsigned int weightIdx = 0;
        double divide_factor = .1;
        w->begin(); 
        while( (weight = w->iterate()) ){
          
          //Divide Factor normalizing weight code
          double normalized_weight = weight->second / minWeight;
          normalized_weight /= divide_factor;
          normalizedIntWeights[weightIdx++] = (int) ceil(normalized_weight); 
          
          //interger normalizing weight code
          // normalizedIntWeights[weightIdx++] = (int) ceil( weight->second / minWeight );
          // divide = false; 

          PCU_Debug_Print("weight %d, normalized to %d from %f\n", weight->first, normalizedIntWeights[(weightIdx-1)], weight->second);
        }

        w->end();            

        const double weightCapacity = (maxW - w->self()); //How much weight can be added on to the current part before it reaches the HeavyImb
        const int knapsackCapacity = floor((weightCapacity/minWeight)/divide_factor); //total weight that can be added to self, normalized to min

        // const int knapsackCapacity = floor(weightCapacity/minWeight);

        PCU_Debug_Print("Weight Capcity = %f \nknapsackCapacity = %d \n", weightCapacity, knapsackCapacity);


        int* value = new int[s->total()];
        std::fill (value, value + s->total(),1);

        knapsack* ks = new knapsack(knapsackCapacity, w->size(), normalizedIntWeights, value);        
        const int solnVal = ks->solve();    
        mergeTargetsResults.reserve(solnVal);
        ks->getSolution(mergeTargetsResults);
        double t1 = MPI_Wtime();
        if (divide) printf("Part %d executed knapsack in %f with knapsack type dividing factor, %f\n", PCU_Comm_Self(), t1 - t0, divide_factor);
        else printf("Part %d executed knapsack in %f with knapsack type interger rounding\n", PCU_Comm_Self(), t1 - t0);
        
        //PCU_Debug_Print("mergetargets start\n");  
        PCU_Debug_Print("mergetargets size = %u\n", mergeTargetsResults.size());

        //Converting mergeTargetsResults to partId's
        vector<int> partIdMergeTargets;
        partIdMergeTargets.reserve(mergeTargetsResults.size());
        
        for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
          partIdMergeTargets.push_back(nborPartIds[mergeTargetsResults[i]]);
          
          PCU_Debug_Print("merge Target result %d -- ", mergeTargetsResults[i]);
          PCU_Debug_Print("merge target part ID %d\n", nborPartIds[mergeTargetsResults[i]]);
        }
        PCU_Debug_Print("\n$$");

        mergeTargetsResults.swap(partIdMergeTargets);



        // for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
        //   PCU_Debug_Print("new merge Target result %d\n", mergeTargetsResults[i]);
        // }
         
        delete [] nborPartIds;
        delete [] value;  
        delete [] normalizedIntWeights;
        delete ks;    
               
      }

      size_t total() {
        return mergeTargetsResults.size();
      }

      const int mergeTargetIndex(int& index){
        return mergeTargetsResults.at(index);
      }

      const int test(){
        return mergeTargetsResults[0];
      }

      
    private:
      MergeTargets();
      vector<int> mergeTargetsResults;
  };

  //Convert part indexes to part ID's in MergeTargets // Done
  //TODO is what i did in mergeTargets with the vector valid?
  //Can i change the return of total to be an int? **Change to size_t
  //how does putting const in front of 158 allow it to be changed.
  //Changed FMDB_mesh_getNeighborPartId to directly push into misLuby part instead of creating a set and then passing it in
  //TODO **write function to pull ith partId from mergeTargetResults**


  apf::Migration* selectMerges(apf::Mesh* m, Sides* s, MergeTargets& tgts) {
    //run MIS and getMergeTargets(...) to determine which 'target'
    // part this part will be merged into then create a Migration 
    // object and set the 'target' part id for each element
    //return the migration object

    
    //Vector to be passed into misLuby (later to be changed to just a single part info since vector is still from multi parts per process)
    vector<misLuby::partInfo> parts;
    parts.reserve(1); //Don't know if this is necessary, thought I'd throw it in anyway, remove if believed to be unnecessary
    
    //Generating misLuby part info for current part 
    misLuby::partInfo part;
    part.id = m->getId(); //same info as FMDB_Part_ID

    //Passing in the adjPartIds (works)
    const Sides::Item* partId; //Kind onf understand the const   
    s->begin();         
    while( (partId = s->iterate()) ) 
      part.adjPartIds.push_back(partId->first);
    s->end();
      /*
      PCU_Debug_Print("adjPartIds:\n");
      vector<int>::iterator itr = part.adjPartIds.begin();
      while(itr != part.adjPartIds.end()){
        PCU_Debug_Print("\t%i\n", *itr++);
        }
      PCU_Debug_Print("adjPartIds end\n");
      */

    int i = 0;
    PCU_Debug_Print("size %f\n", tgts.total());
    // PCU_Debug_Print("test %d\n", tgts.mergeTargetIndex(i));

    // use Sides to get neighbor part ids
    // merging net comes from MergeTargets
    // make variable for random number and set it to 0 for testing purposes
    // we can hardcode the random number in the misLuby::partInfo object for testing purposes
    return new apf::Migration(m);
    //Check to see 
  };

  int splits(Weights* w, double tgtWeight) {
    return static_cast<int>(ceil(w->self()/tgtWeight))-1; 
  }

  int isEmpty(apf::Mesh* m, apf::Migration* plan) {
    return (m->count(m->getDimension()) - plan->count() == 0) ? 1 : 0;
  }

  int totSplits(Weights* w, double tgtWeight) {
    int numSplits = splits(w, tgtWeight);
    PCU_Add_Ints(&numSplits, 1);  // MPI_All_reduce(...,MPI_SUM,...)
    return numSplits;
  }

  int numEmpty(apf::Mesh* m, apf::Migration* plan) {
    int empty = isEmpty(m, plan);
    PCU_Add_Ints(&empty, 1);
    return empty;
  }

  bool canSplit(apf::Mesh* m, Weights* w, apf::Migration* plan, double tgt, int& extra) {
    //PCU_Debug_Print("Empty parts = %f, total Splits = %f\n", numEmpty(m,plan), totSplits(w,tgt));
    extra = numEmpty(m, plan) - totSplits(w, tgt);
    // PCU_Debug_Print("Extra = %d\n", extra);
    if ( extra < 0 ){
      return false;  } 
    else{
      return true;
     }
  }

  double avgWeight(Weights* w) {
    double avg = w->self();
    PCU_Add_Doubles(&avg, 1);
    return avg/PCU_Comm_Peers();
  }

  double maxWeight(Weights* w) {
    double max = w->self();
    PCU_Max_Doubles(&max, 1);
    return max;
  }

  double imbalance(Weights* w) {
    return maxWeight(w)/avgWeight(w);
  }

  double imbalance(apf::Mesh* m, apf::MeshTag* wtag) {
    Sides* s = makeElmBdrySides(m);
    Weights* w = makeEntWeights(m, wtag, s, m->getDimension());
    double imb = imbalance(w);
    delete w; 
    delete s;
    return imb;
  }
  //Should we test the maxWeight first and not do a step first?
  double chi(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w) {
    double testW = maxWeight(w); 
    double step = 0.1 * avgWeight(w); //Not necessarily arbitrary but can be changed.
    bool splits = false;
    int extraEmpties = 0;
    do {
      testW -= step;
      //PCU_Debug_Print("Test Imb = %f\n", testW);
      MergeTargets mergeTgts(s, w, testW);
      // PCU_Debug_Print("Test mergeindex 0, %u",mergeTgts.mergeTargetIndex(i));
      apf::Migration* plan = selectMerges(m, s, mergeTgts); 
      splits = canSplit(m, w, plan, testW, extraEmpties);
      delete plan; // not migrating
   } while ( splits );
   testW += step;
   return testW;
  }

  void split(apf::Mesh* m, Weights* w, double tgt, apf::Migration* plan) {
    const int partId = PCU_Comm_Self();
    int numSplit = splits(w, tgt);
    int empty = isEmpty(m, plan);
    assert(!(numSplit && empty));
    int hl[2] = {numSplit, empty};
    //number the heavies and empties
    PCU_Exscan_Ints(hl, 2);
    //send heavy part ids to brokers
    PCU_Comm_Begin();
    for(int i=0; i<numSplit; i++)
      PCU_COMM_PACK(hl[0]+i, partId);
    PCU_Comm_Send();
    int heavyPartId = 0;
    int count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(heavyPartId);
    }
    assert(count==1);
    //send empty part ids to brokers
    PCU_Comm_Begin();
    if ( empty )
      PCU_COMM_PACK(hl[1], partId); 
    PCU_Comm_Send();
    int emptyPartId = -1;
    count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(emptyPartId);
    }
    assert(count==1);
    //brokers send empty part assignment to heavies
    PCU_Comm_Begin();
    if ( emptyPartId != -1 )
      PCU_COMM_PACK(heavyPartId, emptyPartId); 
    PCU_Comm_Send();
    std::vector<int> tgtEmpties;
    while(PCU_Comm_Listen()) {
      int tgtPartId = 0;
      PCU_COMM_UNPACK(emptyPartId);
      tgtEmpties.push_back(tgtPartId);
    }
    assert( numSplit && tgtEmpties.size() );
    //TODO run async rib 
    //TODO assign rib blocks/sub-parts to tgtEmpties
    //TODO add element empty assignments to plan 
  }

  void hps(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w, double tgt) {
    MergeTargets mergeTargets(s, w, tgt);
    apf::Migration* plan = selectMerges(m, s, mergeTargets);   
    split(m, w, tgt, plan);
    m->migrate(plan);
  }

  class HpsBalancer : public apf::Balancer {
    public:
      HpsBalancer(apf::Mesh* m, int v)
        : mesh(m), verbose(v) 
      {
        (void) verbose; // silence!
      }
      void run(apf::MeshTag* wtag) {
        Sides* sides = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, sides, mesh->getDimension());
        double tgt = chi(mesh, wtag, sides, w);
        PCU_Debug_Print("Final Chi = %f\n",tgt);
        delete sides;
        delete w;
        return; //TODO remove return after testing and put deletes below
        hps(mesh, wtag, sides, w, tgt);

      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        (void) tolerance; // shhh
        double t0 = MPI_Wtime();
        run(weights);
        double elapsed = MPI_Wtime()-t0;
        PCU_Max_Doubles(&elapsed, 1);
        double maxImb = imbalance(mesh, weights);
        if (!PCU_Comm_Self())
          printf("elements balanced to %f in %f seconds\n", maxImb, elapsed);
      }
    private:
      apf::Mesh* mesh;
      int verbose;
  };
}; //end parma namespace

apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, int verbosity) {
  return new parma::HpsBalancer(m, verbosity);
}
