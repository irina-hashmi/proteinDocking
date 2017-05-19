/*

Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: Sep, 2013, Department of Computer Science, George Mason University, VA, USA
Goal: configuration class
Note:
Compile:
Run:
Copyright (c) 2013 Irina Hashmi. All rights reserved.

*/

#ifndef CONFIGURATION_H_INCLUDED
#define CONFIGURATION_H_INCLUDED

#include <vector>
#include <cstring>
#include <memory>

#include "../molecule/atom.h"
#include "../molecule/surface.h"
#include "../molecule/triangle.h"
#include "../utils/transformation.h"
#include "../utils/stringutilities.h"
#include "../utils/fileutilities.h"
#include "../utils/Vector3.h"
#include "../file/pdb.h"
#include "../file/readTriangle.h"
#include "../energy/energy.h"
#include "../math/SampleNeighbor.h"
#include "../math/rmsd1.h"
#include "../math/Random.h"
#include "../math/mymath.h"
#include "../file/readSurface.h"
#include "../utils/MyTimer.h"

#include "../learning/Classifier1.h"
#include "../learning/Classifier2.h"
#include "../learning/Classifier3.h"
#include "../learning/Classifier4.h"
#include "../learning/Classifier5.h"
#include "../learning/Classifier6.h"
#include "../learning/Classifier7.h"
#include "../learning/Classifier8.h"
#include "../learning/Classifier9.h"
#include "../learning/Classifier10.h"


namespace ihashmi{

  class Configuration{
    public:


      /*constructor & destructor*/
      Configuration();
      //Configuration(Configuration &c);
      ~Configuration();
      void PrintAConfiguration();
      void WriteConfigurationInPDBFormat();
      void WriteConfigurationInPDBFormatInLogFile(std::string fileName,
                                                   int generationNo,
                                                   int populationNo);

      void WriteConfigurationInLogFile(ofstream& logFile,
                                       int generationNo,
                                       int individualNo);

      void WriteConfigurationInCoordinateFile(ofstream& coordinateFile,
                                              int generationNo,
                                              int individualNo);

      void InitializeConfiguration(std::string pdbId_,
                                    std::string referenceChain_,
                                    std::string moveChain_,
                                    std::string chains_);

      // update the triangle, transformation info
      void UpdateAfterDocking(std::string pdbId_, std::string allChains_,
                              std::string moveChain_, std::string refChain_,
                              unsigned long moveTriId_, unsigned long refTriId_,
                              Vector3 translationVector_, real* rotationMatrix_,
                              unsigned long refAtomStartIndex_,
                              unsigned long refAtomEndIndex_,
                              unsigned long moveAtomStartIndex_,
                              unsigned long moveAtomEndIndex_,
                              real referenceConfiguraionInteractionEnergy_);

      void PrepareForDocking(); // calculate the ms/cp/triangle

      //pdbId
      void SetPdbId(std::string pdbId_);
      std::string GetPdbId();

      /* ****************** native, ms = cp and triangle****************** */
      void CalculateConnollyRepresentation();
      void CalculateShuoRepresentation();
      void CalculateTriangleRepresentation();

      void SetCurrentConfigurationInfo(std::vector <Atom>& v);
      void GetCurrentConfigurationInfo(std::vector <Atom>& v);
      void PrintCurrentConfigurationInfo(std::string fileName);
      unsigned long GetCurrentConfigurationSize();

      void SetConnollyPointInfo(std::vector <Surface>& s);
      void GetConnollyPointInfo(std::vector <Surface>& s);
      void PrintConnollyPointInfo(std::string fileName);
      unsigned long GetConnollyPointInfoSize();


      void SetShuoPointInfo(std::vector <Surface>& s);
      void GetShuoPointInfo(std::vector <Surface>& s);
      void PrintShuoPointInfo(std::string fileName);
      unsigned long GetShuoPointInfoSize();


      void SetTriangleInfo(std::vector <Triangle>& t);
      void GetTriangleInfo(std::vector <Triangle>& t);
      void PrintTriangleInfoVector(std::string fileName);
      int GetTriangleInfoVectorSize();

      /* ****************** dock and transformation ****************** */
      bool CopyDock(Configuration& c1,
                           Configuration& c2,
                           std::vector <Configuration> &v);

      unsigned long CopyTransformConfiguration(Configuration& c2,
                                           real* finalRotationMatrix,
                                           Vector3& finalTranslationVector3);

      //dock is for extension
      bool Dock(Configuration& c,
                           Configuration& dockedConfig,
                           std::vector <Configuration> &v);

      void DockAroundNeighborhood(std::vector <Triangle> &t1,
                                             std::vector <Triangle> &t2,
                                             int neighborTriangleId1,
                                             int neighborTriangleId2,
                                             Configuration& c1,
                                            Configuration& c2,
                                            std::vector <Configuration> &v);

      void RandomDock(Configuration& c2,
                      Configuration& dockedConfig,
                      std::vector <Configuration>& v);


      //pass a whole configuration but transform from start index to end index
      unsigned long TransformConfigurationWithIndex(real *finalRotationMatrix,
                  Vector3& finalTranslationVector3, Configuration& dockedConfig,
                  unsigned long start, unsigned long end);

      unsigned long TransformConfiguration(Configuration& c,
                  real *finalRotationMatrix, Vector3& finalTranslationVector3,
                  Configuration& dockedConfig);
      //c to be randomize and dock again used for modification
      //bool RandomDock(Configuration *dockedConfig, std::vector <Configuration*>& v);

      void GetRotationMatrix(real *resultMatrix);
      Vector3 GetTranslationVector();
      void SetRotationMatrix(real *rotationMatrix_);
      void SetTranslationVector(Vector3& v);

      /* ****************** DT ****************** */
      int IfTrueInteractionOld();
      int IfTrueInteractionNew();


      /* ****************** chains ****************** */
      void SetAllChains(std::string chains);
      std::string GetAllChains();
      void SetReferenceChainId(std::string chain);
      void SetMoveChainId(std::string chain);
      std::string GetReferenceChainId();
      std::string GetMoveChainId();
      std::string GetRemainingUnit();
      void SetRemainingUnit();
      std::string GetInvolvedChainId();

      /* ************* parent id in config vector ****************** */
      //void SetReferenceConfigurationIndex(unsigned long refConfigIndex_)     ;
      //unsigned long GetReferenceConfigurationIndex();

      /* ************* interface and triangle information ****************** */
      void SetInterfaceAtoms();
      void GetInterfaceAtoms(std::vector <int>& interface);
      unsigned long GetReferenceTriangleId();
      unsigned long GetMoveTriangleId();
      void SetReferenceTriangleId(unsigned long id);
      void SetMoveTriangleId(unsigned long id);
      unsigned long referenceAtomStartIndex;
      unsigned long referenceAtomEndIndex;
      unsigned long moveAtomStartIndex;
      unsigned long moveAtomEndIndex;
      void DetermineInterfaceAtoms();
      std::vector <int> referenceAtoms;
      std::vector <int> moveAtoms;
      friend class EA;
      friend class Generation;

      /* ***************** energy and score, COM *********************** */
      real CalculateTotalInteractionEnergy();
      void SetTotalInteractionEnergy(real totalInteractionEnergy_);
      const real GetTotalInteractionEnergy() const;
      real GetReferenceConfigurationInteractionEnergy();
      void SetReferenceConfigurationInteractionEnergy(real
                          referenceConfigurationInteractionEnergy_);

      void SetScaledInteractionEnergy(real val_);
      real GetScaledInteractionEnergy(void);

      real CalculatelowestlRMSD(std::vector <Configuration> &v);
      real MeasureSimilarity(Configuration &c);
      void SetlowestlRMSD(real lowestlRMSD_);
      real GetlowestlRMSD();
      int generationNo, individualNo;
      Vector3 CalculateCenterOfMass(std::vector <Atom>&a);
      int GetDTScore();
      void SetDTScore(int val_);


      //update configuration
      //void UpdateConfiguration();

      //misc util files
      //copy config/atom info from a to this starting from start to end
      void CopyConfigurationInfo(std::vector <Atom>&src, unsigned long start,
                            unsigned long end);

      unsigned long GetParentId();
      void SetParentId(unsigned long parentId_);




      //clustering variables
      int GetMembershipStatus(); //0 = not selected, 1 = leader, 2 = memeber
      void SetMembershipStatus(int status_); //0 = not selected, 1 = leader, 2 = memeber

      const unsigned long GetClusterSize() const; //applicable for leaders only
      void SetClusterSize(int size_);

      void IncreaseClusterSize(void);


      /* ********************* d2d stuff *********************** */
      bool d2dDock(Configuration& c,
                    Configuration& dockedConfig,
                    std::vector <Configuration> &v, int totalConfigs);

      void CalculateSpecificShuoRepresentation(); //only focus on TM regions
      void PrepareForDockingInSpecificRegion();
      real GetTotalTimingInfo();
      void SetTotalTimingInfo(real time_);

      real GetEnergyTimingInfo();
      void SetEnergyTimingInfo(real time_);







    private:

      std::string pdbId;
      //real foldXBindingEnergy;
      string remainingUnit;
      real rotationMatrix[9];
      Vector3 translationVector;
      string allChains;
      string referenceChainId;
      string moveChainId;
      unsigned long referenceTriangleId; //reference = receptor
      unsigned long moveTriangleId;      //move = ligand
      real lowestlRMSD;
      real totalInteractionEnergy;
      real scaledInteractionEnergy;
      unsigned long parentId; //parent id is the reference config id in configuration vector
      //Vector3 centerOfMass;

      real referenceConfigurationInteractionEnergy;
      std::vector <Atom> currentConfigurationInfo; //current points
      std::vector <Surface> connollyPointInfoVector; //connolly/ ms point
      std::vector <Surface> shuoPointInfoVector; //critical point
      std::vector <Triangle> triangleInfoVector; //triangle


      //clustering variables
      int membershipStatus; //0 = not selected, 1 = leader, 2 = memeber
      unsigned long clusterSize; //applicable for leaders only
      int whichCluster; //which cluster you are in
      real avgIntraClusterDist; //later
      int DTscore;
      real timingInfo;
      real energyTimingInfo;
      //unsigned long referenceConfigurationIndex;

  }; //class
} //namespace


#endif // CONFIGURATION_H_INCLUDED
