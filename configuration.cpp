/*

Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: Sep, 2013, Department of Computer Science, George Mason University, VA, USA
Goal: main.cpp
Note:
Compile:
Run:
Copyright (c) 2013 Irina Hashmi. All rights reserved.

*/

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "configuration.h"
#include "../utils/matrix/matrix.h"
#include "../utils/definitions.h"
#include <memory>

//using namespace std;
namespace ihashmi{
  //default constructor
  Configuration::Configuration(){
    //initialize

  }
  //default destructor
  Configuration::~Configuration(){
    //delete a configuration
    /*for(int i = 0; i < GetCurrentConfigurationSize(); i++){
      delete currentConfigurationInfo[i];

    }
    for(int i = 0; i < GetConnollyPointInfoSize(); i++){
      delete connollyPointInfoVector[i];

    }

    for(int i = 0; i < GetShuoPointInfoSize(); i++){
      delete shuoPointInfoVector[i];

    }

    for(int i = 0; i < GetTriangleInfoVectorSize(); i++){
      delete triangleInfoVector[i];

    }*/
    //destroy();

  }

  void Configuration::WriteConfigurationInPDBFormat(){

    ofstream oFile;
    std::string fileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";
    oFile.open(StringToChar(fileName));
    //cout << "filename in WriteConfig: " << fileName << endl;
    string spacer(" ");
    int dash;

    for(int j = 0 ; j < GetCurrentConfigurationSize();j++){
      stringstream atom_stream;
      if (currentConfigurationInfo[j].GetAtomName().size() > 3){
          dash = 0;
      }
      else{
        dash = 1;
      }
      atom_stream << setiosflags(ios::left) << (spacer.c_str()+1)-dash
                  << setw(4-dash) << currentConfigurationInfo[j].GetAtomName()
                  << " ";
      oFile << setiosflags(ios:: right)
            << setiosflags(ios::fixed) << setprecision(THREE);
      oFile << "ATOM  " << setw(FIVE)
            << currentConfigurationInfo[j].GetAtomSerialNumber() <<  " "
            << atom_stream.str() << currentConfigurationInfo[j].GetResName();
      oFile << " " << currentConfigurationInfo[j].GetChainID() << " "
            << setw(THREE)<< currentConfigurationInfo[j].GetResSeq() << "    ";
      oFile << setw(EIGHT) << currentConfigurationInfo[j].GetX()
            << setw(EIGHT) << currentConfigurationInfo[j].GetY()
            << setw(EIGHT) << currentConfigurationInfo[j].GetZ();

      oFile << "  " << "1.00" << endl;
    }
    oFile << "TER" << endl;
    oFile.close();
  }
  //Generation # Individual # Chains Energy RMSD
  //coordinate information in the following format
  //x1 y1, z1 \n
  //x2 y2, z2 \n
  // ........
  //xn yn, zn \n
  void Configuration::WriteConfigurationInCoordinateFile
  (ofstream& coordinateFile, int genNo,int indvNo){

    coordinateFile
          << "Generation " << generationNo << " Individual " << indvNo
          << " Chains " << GetInvolvedChainId()
          << " Size " << GetCurrentConfigurationSize()
          << " Energy " << GetTotalInteractionEnergy()
          << " RMSD " << GetlowestlRMSD() << endl;

    for(int j = 0 ; j < GetCurrentConfigurationSize();j++){
      coordinateFile << "\t" << currentConfigurationInfo[j].GetX()
                     << "\t" << currentConfigurationInfo[j].GetY()
                     << "\t" << currentConfigurationInfo[j].GetZ();
    }
    coordinateFile << "\nTER\n" << endl;

  }

  //Gen # Indv # Chains Size lRMSD energy
  void Configuration::WriteConfigurationInLogFile(ofstream& logFile,int genNo,int indvNo){
    logFile << genNo << "\t" << indvNo << "\t" << GetInvolvedChainId() << "\t";
    logFile << GetCurrentConfigurationSize() << "\t" << GetlowestlRMSD() << "\t";
    logFile << GetTotalInteractionEnergy() << endl;

  }


  void Configuration::InitializeConfiguration(std::string pdbId_,
                                              std::string refChain_,
                                              std::string moveChain_,
                                              std::string chains_){

    SetPdbId(pdbId_);
    SetReferenceChainId(refChain_);
    SetMoveChainId(moveChain_);
    SetAllChains(chains_);
    SetMembershipStatus(0);
    SetClusterSize(0);
  }

  void Configuration::PrepareForDocking(){
    WriteConfigurationInPDBFormat();
    CalculateConnollyRepresentation();
    CalculateShuoRepresentation();
    CalculateTriangleRepresentation();
  }



  void Configuration::UpdateAfterDocking(std::string pdbId_,
                                 std::string allChains_,
                                 std::string moveChain_,
                                 std::string refChain_,
                                 unsigned long moveTriId_,
                                 unsigned long refTriId_,
                                 Vector3 translationVector_,
                                 real* rotationMatrix_,
                                 unsigned long refAtomStartIndex_,
                                 unsigned long refAtomEndIndex_,
                                 unsigned long moveAtomStartIndex_,
                                 unsigned long moveAtomEndIndex_,
                                 //unsigned long referenceConfigurationIndex_,
                                 real referenceConfiguraionInteractionEnergy_
                                 /*unsigned long parentId_
                                 real totalInteractionEnergy_*/){
    SetPdbId(pdbId_);
    SetAllChains(allChains_);
    SetMoveChainId(moveChain_);
    SetReferenceChainId(refChain_);

    SetMoveTriangleId(moveTriId_);
    SetReferenceTriangleId(refTriId_);
    SetTranslationVector(translationVector_);
    SetRotationMatrix(rotationMatrix_);
    //SetTotalInteractionEnergy(totalInteractionEnergy_);
    //SetReferenceConfigurationIndex(referenceConfigurationIndex_);
    SetReferenceConfigurationInteractionEnergy(referenceConfiguraionInteractionEnergy_);
    referenceAtomStartIndex = refAtomStartIndex_;
    referenceAtomEndIndex = refAtomEndIndex_;
    moveAtomStartIndex = moveAtomStartIndex_;
    moveAtomEndIndex = moveAtomEndIndex_;
    //SetParentId(parentId_);

    //SetCurrentConfigurationInfo(currentConfigurationInfo_);
    //SetCurrentConfigurationInfo(currentConfigVect_);
  }

  void Configuration::SetPdbId(std::string pdbId_){
    pdbId = pdbId_;
  }

  std::string Configuration::GetPdbId(){
    return pdbId;
  }

  void Configuration::CalculateConnollyRepresentation(){

    std::string pdbFileName = this->GetPdbId() + "_" +
      this->GetInvolvedChainId() + ".pdb";
    std::string str = "./runMSPoints_edit " + pdbFileName + " > msout";
    const char* cmd = StringToChar(str);
    system(cmd);

  }

  void Configuration::CalculateShuoRepresentation(){

    std::string pdbFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";
    std::string msFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb.ms";
    std::string shuoFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb.shuo";
    std::string str = "./shuo_surface.Linux "
          + pdbFileName + " " + msFileName + " " + shuoFileName  + " > shuoout";
    const char* cmd = StringToChar(str);
    system(cmd);

    //keep the pdb file for the energy evaluation
    //remove MS file
    RemoveFile(msFileName);

  }

  void Configuration::CalculateTriangleRepresentation(){

    std::string shuoFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb.shuo";
    std::string triangleFileName = GetPdbId() + "_" + GetInvolvedChainId()
                                         + ".unique.triangle.score.out";
    //generates the .unique.triangle.score.out file
    std::string str = "python generateTriangle.py " + GetPdbId() + ' '
                                                    + GetInvolvedChainId();
    const char* cmd = StringToChar(str);
    system(cmd);
    std::vector<Triangle> t;
    unsigned long totalTriangle = ReadTriangle(triangleFileName, t);
    SetTriangleInfo(t);
    //remove .shuo.pdb and .unique.triangle.score.out
    //RemoveFile(shuoFileName);
    //RemoveFile(triangleFileName);
    t.clear();
  }


  /* ********************* connolly/shuo/atoms ******************* */
  //set the current atomic info
  //initialize with native info
  void Configuration::SetCurrentConfigurationInfo(std::vector <Atom>& v){
    currentConfigurationInfo = v;

  }
  //get the current atomic info
  //v is the holder for getting the current configuration
  void Configuration::GetCurrentConfigurationInfo(std::vector <Atom>& v){
    v = currentConfigurationInfo;

  }

  unsigned long Configuration::GetCurrentConfigurationSize(){
    return currentConfigurationInfo.size();

  }

  void Configuration::SetConnollyPointInfo(std::vector <Surface>& s){
    connollyPointInfoVector = s;
  }
  void Configuration::GetConnollyPointInfo(std::vector <Surface>& s){
    s = connollyPointInfoVector;

  }

  unsigned long Configuration::GetConnollyPointInfoSize(){
    return connollyPointInfoVector.size();

  }

  void Configuration::SetShuoPointInfo(std::vector <Surface>& s){
    this->shuoPointInfoVector = s;

  }
  void Configuration::GetShuoPointInfo(vector <Surface>& s){
    s = this->shuoPointInfoVector;


  }
  unsigned long Configuration::GetShuoPointInfoSize(){
    return shuoPointInfoVector.size();

  }

  /* ********************* triangle ******************* */

  void Configuration::SetTriangleInfo(std::vector <Triangle>& t){
    this->triangleInfoVector = t;
  }

  void Configuration::GetTriangleInfo(std::vector <Triangle>& t){
    t = this->triangleInfoVector;

  }

  int Configuration::GetTriangleInfoVectorSize(){
    return this->triangleInfoVector.size();

  }

  /* ********************* Rotation/Translation component ******************* */

  void Configuration::GetRotationMatrix(real *rotationMatrix_){
      CopyMatrix(rotationMatrix_, rotationMatrix);

  }

  Vector3 Configuration::GetTranslationVector(){
    return translationVector;

  }

  void Configuration::SetRotationMatrix(real *rotationMatrix_){
    CopyMatrix(rotationMatrix, rotationMatrix_);

  }

  void Configuration::SetTranslationVector(Vector3& v){
      this->translationVector.update(v);

  }

  /* *************** Interface, triangle and energy ********************** */
  /*void Configuration::DetermineInterfaceAtoms(){
    for(int i = this->referenceAtomStartIndex;
      i < this->referenceAtomEndIndex; i++ ){
      for(int j = this->moveAtomStartIndex;
        j < this->moveAtomEndIndex; j++)

        if (this->currentConfigurationInfo[i]->
            EucledianDistanceSqr(this->currentConfigurationInfo[j])
                                                        <= INTERFACETHRESHOLD){
            this->referenceAtoms.push_back(
                    this->currentConfigurationInfo[i]->GetAtomSerialNumber());
            this->moveAtoms.push_back(
                    this->currentConfigurationInfo[j]->GetAtomSerialNumber());
        }

    }

  }*/

  void Configuration::SetReferenceTriangleId(unsigned long id){
      referenceTriangleId = id;
  }

  void Configuration::SetMoveTriangleId(unsigned long id){
    moveTriangleId = id;

  }

  unsigned long Configuration::GetReferenceTriangleId(){
        return referenceTriangleId;
  }

  unsigned long Configuration::GetMoveTriangleId(){
      return moveTriangleId;

  }

  /*real Configuration::CalculateTotalInteractionEnergy(bool flag, real parentTotalEnergy){
    Energy *e = new Energy();
    std::string pdbFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";
    real totalEnergy=0;
    real currInteractionEnergy = e->GetFoldXInteractionEnergy(pdbFileName,
                                                   GetReferenceChainId(),
                                                   GetMoveChainId());
    if(GetInvolvedChainId().size()==2){//
      totalEnergy = currInteractionEnergy;
      cout << "1. in calculate energy " << totalEnergy << endl;
    }
    else{
      if(GetInvolvedChainId().size() > 2 and flag == true ){//extension flag
        totalEnergy = GetReferenceConfigurationInteractionEnergy() +
                                currInteractionEnergy;
        cout << "2. in calculate energy " << totalEnergy << endl;
      }
      else{//total chain > 2 and modify
        totalEnergy = GetReferenceConfigurationInteractionEnergy() +
                      currInteractionEnergy - (parentTotalEnergy -        GetReferenceConfigurationInteractionEnergy());
        cout << "3. in calculate energy " << totalEnergy << endl;
      }

    }
    delete e;
    return totalEnergy;
  }*/

  void Configuration::SetScaledInteractionEnergy(real val_){
    scaledInteractionEnergy = val_;
  }


  real Configuration::GetScaledInteractionEnergy(void){
    return scaledInteractionEnergy;


  }

  real Configuration::CalculateTotalInteractionEnergy(){
    Energy *e = new Energy();
    std::string pdbFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";
    real totalEnergy=0;
    real currInteractionEnergy = e->GetFoldXInteractionEnergy(pdbFileName,
                                                   GetReferenceChainId(),
                                                   GetMoveChainId());

    std::string c = this->GetInvolvedChainId();
    int n = c.size();
    //calculate nC2
    if(FITNESS==1){
      totalEnergy = currInteractionEnergy/Choose(n,2);
    }
    else{
      totalEnergy = currInteractionEnergy;
    }
    delete e;
    return totalEnergy;
  }


  void Configuration::SetTotalInteractionEnergy(real interactionEnergy_){
    totalInteractionEnergy = interactionEnergy_;

  }

  const real Configuration::GetTotalInteractionEnergy() const{
    return totalInteractionEnergy;

  }

  real Configuration::GetReferenceConfigurationInteractionEnergy(){
    return referenceConfigurationInteractionEnergy;
  }

  void Configuration::SetReferenceConfigurationInteractionEnergy(real referenceConfigurationInteractionEnergy_){
    referenceConfigurationInteractionEnergy =
            referenceConfigurationInteractionEnergy_;

  }

  //similarity measurement with this to passed configurations
  //SIMILARITYMEASURE == RMSD
  //SIMILARITYMEASURE == INTERFACE (to be implemented)
  real Configuration::MeasureSimilarity(Configuration &c){
    if (SIMILARITYMEASURE){
      //prepare the native array
      double x[GetCurrentConfigurationSize()*3]; // native
      double y[GetCurrentConfigurationSize()*3]; //current

      unsigned long k = 0;

      for(unsigned long i = 0; i < GetCurrentConfigurationSize(); i++){
        y[k] = c.currentConfigurationInfo[i].GetX();
        y[k+1] = c.currentConfigurationInfo[i].GetY();
        y[k+2] = c.currentConfigurationInfo[i].GetZ();
        k+=3;
      }

      k = 0;
      for(unsigned long i = 0; i < GetCurrentConfigurationSize(); i++){
        x[k] = currentConfigurationInfo[i].GetX();
        x[k+1] = currentConfigurationInfo[i].GetY();
        x[k+2] = currentConfigurationInfo[i].GetZ();
        k+=3;
      }
      real lRMSD = rmsd(GetCurrentConfigurationSize()*3,x,y);
      return lRMSD;

    }//RMSD
    else{
      cout << "To be implemented:" << endl;
      cout << "Provision for contact interface similarity measure: dockrank" << endl;
      cout << "or using pareto optimization for similarity. " << endl;
      cout << "or lrmsd and interface. " << endl;

    }
  }


  void Configuration::SetlowestlRMSD(real lowestlRMSD_){
    lowestlRMSD = lowestlRMSD_;
  }


  real Configuration::GetlowestlRMSD(){
    return lowestlRMSD;
  }
  //this is a slopy code I need to chage the original rmsd code
  real Configuration::CalculatelowestlRMSD(std::vector <Configuration> &v){
    std::string s = GetInvolvedChainId();
    std::vector<unsigned long> chainVector;
    //cout << "inside rmsd calculation:chains are: " << s << endl;
    for(int i = 0 ; i < s.size(); i++){
      chainVector.push_back(MapCharToDigit(s[i]));
      //cout << "chainVector " << chainVector[i] << endl;
    }
    //prepare the native array
    double x[GetCurrentConfigurationSize()*3]; // native
    double y[GetCurrentConfigurationSize()*3]; //current

    //cout << "size of x and y " << GetCurrentConfigurationSize()*3 << endl;

    unsigned long k = 0;
    //native configuration into array x for rmsd
    //cout << "native points\n";
    for(int l = 0; l < chainVector.size(); l++){
      for(int i = 0 ; i < v[chainVector[l]].GetCurrentConfigurationSize(); i++){
        x[k] = v[chainVector[l]].currentConfigurationInfo[i].GetX();
        x[k+1] = v[chainVector[l]].currentConfigurationInfo[i].GetY();
        x[k+2] = v[chainVector[l]].currentConfigurationInfo[i].GetZ();
        //cout << x[k] << "\t" << x[k+1] << "\t" << x[k+2] << "\t";
        k+=3;
      }//i

    }//l
    //cout << "End Native\n";

    k=0;
    //current configuration into array y for rmsd
    for(unsigned long i = 0; i < GetCurrentConfigurationSize(); i++){
      y[k] = currentConfigurationInfo[i].GetX();
      y[k+1] = currentConfigurationInfo[i].GetY();
      y[k+2] = currentConfigurationInfo[i].GetZ();
      //cout << y[k] << "\t" << y[k+1] << "\t" << y[k+2] << "\t";
      k+=3;
    }
    //cout << "End of computed structure\n";

    //cout  << "\n\nsize of x and y " << (sizeof x / sizeof x[0]) << "\t"
      //    << (sizeof y / sizeof y[0]) << endl;

    real lRMSD = rmsd(GetCurrentConfigurationSize()*3,x,y);
    SetlowestlRMSD(lRMSD);
    //cout << "end lrmsd calculateion: rmsd " << lRMSD << endl;
    return lRMSD;
  }


  /* ********************** Chains ****************************** */
  std::string Configuration::GetInvolvedChainId(){
    if ( GetMoveChainId().compare("NULL") != 0)
      return GetReferenceChainId() + GetMoveChainId();
    else{
      return GetReferenceChainId();
    }
  }

  void Configuration::SetAllChains(std::string chain){
    allChains = chain;

  }

  std::string Configuration::GetAllChains(){
    return allChains;
  }

  void Configuration::SetReferenceChainId(std::string chain){
    referenceChainId = chain;

  }
  void Configuration::SetMoveChainId(std::string chain){
    moveChainId = chain;

  }
  std::string Configuration::GetReferenceChainId(){
    return referenceChainId;

  }
  std::string Configuration::GetMoveChainId(){
    return moveChainId;

  }
  //remaining units
  std::string Configuration::GetRemainingUnit(){
    std::string involvedChain = GetReferenceChainId() + GetMoveChainId();
    return difference(GetAllChains(), involvedChain);
  }


  /* ******************* Dock and Transformation ********************** */
  //transform part of a configurations
  //specially written for d2d
  //where we transform only from the start to end
  unsigned long Configuration::TransformConfigurationWithIndex
                  (real *finalRotationMatrix,
                  Vector3& finalTranslationVector3,
                  Configuration& dockedConfig,
                  unsigned long start,
                  unsigned long end){

    for(unsigned long k = start; k <= end; k++){

        Vector3 movePrevPoint(this->currentConfigurationInfo[k].GetX(),
                              this->currentConfigurationInfo[k].GetY(),
                              this->currentConfigurationInfo[k].GetZ());
        Vector3 rotateMovePrevPoint =
                          MultMatrixVector(finalRotationMatrix, movePrevPoint);
        Vector3 transformedMovePrevPoint =
                            rotateMovePrevPoint + finalTranslationVector3;
        Atom tempAtom
               (this->currentConfigurationInfo[k].GetAtomSerialNumber(),
                this->currentConfigurationInfo[k].GetAtomName(),
                this->currentConfigurationInfo[k].GetAltLoc(),
                this->currentConfigurationInfo[k].GetResName(),
                this->currentConfigurationInfo[k].GetChainID(),
                this->currentConfigurationInfo[k].GetResSeq(),
                transformedMovePrevPoint[0],
                transformedMovePrevPoint[1],
                transformedMovePrevPoint[2],
                this->currentConfigurationInfo[k].GetOccupancy(),
                this->currentConfigurationInfo[k].GetTempFactor(),
                this->currentConfigurationInfo[k].GetElement(),
                false);
        //cout << this->currentConfigurationInfo[k]->GetAtomSerialNumber() < "\t";
        dockedConfig.currentConfigurationInfo.push_back(tempAtom);
    }
    return dockedConfig.GetCurrentConfigurationSize();
  }

  //transform c and push back to dockedConfig used by Dock method
  unsigned long Configuration::TransformConfiguration(Configuration& c,
                  real *finalRotationMatrix, Vector3& finalTranslationVector3,
                  Configuration& dockedConfig){
    for(int k = 0; k < c.GetCurrentConfigurationSize(); k++){

        Vector3 movePrevPoint(c.currentConfigurationInfo[k].GetX(),
                              c.currentConfigurationInfo[k].GetY(),
                              c.currentConfigurationInfo[k].GetZ());
        Vector3 rotateMovePrevPoint =
                          MultMatrixVector(finalRotationMatrix, movePrevPoint);
        Vector3 transformedMovePrevPoint =
                            rotateMovePrevPoint + finalTranslationVector3;
        Atom tempAtom
               (c.currentConfigurationInfo[k].GetAtomSerialNumber(),
                c.currentConfigurationInfo[k].GetAtomName(),
                c.currentConfigurationInfo[k].GetAltLoc(),
                c.currentConfigurationInfo[k].GetResName(),
                c.currentConfigurationInfo[k].GetChainID(),
                c.currentConfigurationInfo[k].GetResSeq(),
                transformedMovePrevPoint[0],
                transformedMovePrevPoint[1],
                transformedMovePrevPoint[2],
                c.currentConfigurationInfo[k].GetOccupancy(),
                c.currentConfigurationInfo[k].GetTempFactor(),
                c.currentConfigurationInfo[k].GetElement(),
                false);
        //cout << this->currentConfigurationInfo[k].GetAtomSerialNumber() < "\t";
        dockedConfig.currentConfigurationInfo.push_back(tempAtom);
        //delete tempAtom;

    }
    //cout << endl;
    //dockedConfig.PrintAConfiguration();
    return dockedConfig.GetCurrentConfigurationSize();
  }

  unsigned long Configuration::CopyTransformConfiguration(
                                           Configuration& c2,
                                           real* finalRotationMatrix,
                                           Vector3& finalTranslationVector3){

     for(int k = 0; k < c2.GetCurrentConfigurationSize(); k++){
        Vector3 movePrevPoint(c2.currentConfigurationInfo[k].GetX(),
                              c2.currentConfigurationInfo[k].GetY(),
                              c2.currentConfigurationInfo[k].GetZ());
        Vector3 rotateMovePrevPoint =
                          MultMatrixVector(finalRotationMatrix, movePrevPoint);
        Vector3 transformedMovePrevPoint =
                            rotateMovePrevPoint + finalTranslationVector3;
        Atom tempAtom
               (c2.currentConfigurationInfo[k].GetAtomSerialNumber(),
                c2.currentConfigurationInfo[k].GetAtomName(),
                c2.currentConfigurationInfo[k].GetAltLoc(),
                c2.currentConfigurationInfo[k].GetResName(),
                c2.currentConfigurationInfo[k].GetChainID(),
                c2.currentConfigurationInfo[k].GetResSeq(),
                transformedMovePrevPoint[0],
                transformedMovePrevPoint[1],
                transformedMovePrevPoint[2],
                c2.currentConfigurationInfo[k].GetOccupancy(),
                c2.currentConfigurationInfo[k].GetTempFactor(),
                c2.currentConfigurationInfo[k].GetElement(),
                false);
        //cout << this->currentConfigurationInfo[k].GetAtomSerialNumber() < "\t";
        this->currentConfigurationInfo.push_back(tempAtom);

    }
    //cout << endl;
    return this->GetCurrentConfigurationSize();
  }




  //will only be called from expansion
  //c1 and c2 are reference and moving configuration
  //after docking c1 and c2 put it in the dockedConfig
  //this reference referes to dockedConfig
  bool Configuration::CopyDock(Configuration& c1,
                           Configuration& c2,
                           std::vector <Configuration> &v){

    //randomly sample two triangles from the triangle vector
    /*cout << "Inside Copy Dock: ";
    cout << "before docking: ";
    cout << "curr reference and move unit triangle size: "<<endl;
    cout << c1.GetCurrentConfigurationSize() << endl;
    cout << c2.GetTriangleInfoVectorSize() << endl;
    cout << this->GetCurrentConfigurationSize() << endl;*/

    int referenceTriangleId = IntUrandom(this->GetTriangleInfoVectorSize());
    int moveTriangleId = IntUrandom(c2.GetTriangleInfoVectorSize());

    //find two geometrically fit triangle
    //if fit then transform
    int try1 = TRIANGLETRYCOUNT;
    int try2 = TRIANGLETRYCOUNT;

    bool found1 = false;
    bool found2 = false;
    bool found = false;

    while(try1 > 0  and found1 != true){

      int neighborTriangleId1 = findNeighborTriangle(referenceTriangleId, this->triangleInfoVector);
      if(neighborTriangleId1 < 0){
        try1--;
        continue;

      }
      while(try2 > 0 and  found2 != true){
        int neighborTriangleId2 = findNeighborTriangle(moveTriangleId, c2.triangleInfoVector);
        if(neighborTriangleId2 < 0){
          try2--;
          continue;

        }
        //if geometric shape complement and other geometric constraints meet
        if(triangleInfoVector[neighborTriangleId1].CheckShapeComplementarity
          (c2.triangleInfoVector[neighborTriangleId2])
          &&  (triangleInfoVector[neighborTriangleId1].PairwiseTriangle
            (c2.triangleInfoVector[neighborTriangleId2]))){
            //if shape and geometrically complementary
            //transform move triangle on top of reference
            //obtain final rotation and translation
            real referenceRotationMatrix[9], moveRotationMatrix[9];
            Vector3 referenceTranslationVector3, moveTranslationVector3;

            //relative rotation and translation of reference frame wrt global
            this->triangleInfoVector[neighborTriangleId1].TransformTriangle
                  (referenceRotationMatrix, referenceTranslationVector3);
            //relative rotation and translation of move frame wrt global
            c2.triangleInfoVector[neighborTriangleId2].TransformTriangle
                  (moveRotationMatrix, moveTranslationVector3);
            //update triangle info
            this->SetReferenceTriangleId(neighborTriangleId1);
            c2.SetMoveTriangleId(neighborTriangleId2);

            //calculate the final translation and rotation
            //final rotation
            //if two traingles are P(reference) and Q(move) wrt global
            //final rot from Q to P . R_PQ = R_WP.transpose() * R_WQ
            real referenceRotationMatrixTranspose[9];
            real finalRotationMatrix[9];
            //1. R_WP.transpose()
            TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
            //2. R_WP.transpose() * R_WQ
            MultMatrix(finalRotationMatrix,
                       referenceRotationMatrixTranspose,
                       moveRotationMatrix);

            //final translation
            //T = R_WP * A - R_WQ*D
            //T_PQ = R_WP.transpose * T
            Vector3 distanceOrigReferenceMove =
                          moveTranslationVector3 - referenceTranslationVector3;
            Vector3 finalTranslationVector3 =
            MultMatrixVector(referenceRotationMatrixTranspose,
                             distanceOrigReferenceMove);

            //transform move configuration on top of refernce configuration
            //using reference and move frame rotation
            //applying on move frame
            //return the whole docked configuration in dockedConfig object

            //transform c2 and push_back to dockedConfig
            unsigned long l = this->CopyTransformConfiguration
                  (c2, finalRotationMatrix,
                   finalTranslationVector3);

            this->UpdateAfterDocking
                  (c1.GetPdbId(), c2.GetAllChains(),
                  c2.GetInvolvedChainId(), c1.GetInvolvedChainId(),
                  neighborTriangleId2, neighborTriangleId1,
                  finalTranslationVector3, finalRotationMatrix,
                  0, c1.GetCurrentConfigurationSize()-1,
                  c1.GetCurrentConfigurationSize(),
                  this->GetCurrentConfigurationSize()-1,
                  c1.GetReferenceConfigurationInteractionEnergy());
                  /*dockedConfig.CalculateTotalInteractionEnergy(true,
                    this->GetTotalInteractionEnergy())*/
            /*cout << "NeighborTriangleIds: "
                  << neighborTriangleId1 << "\t" << neighborTriangleId2 << endl;
            cout << "Updated value after docking: ";
            cout << this->GetInvolvedChainId() << "\t"
                  << this->GetReferenceChainId() << "\t"
                  << this->GetMoveChainId() << "\t"
                  //<< this->GetParentId() << "\t"
                  << this->GetReferenceTriangleId() << "\t"
                  << this->GetMoveTriangleId() << "\t"
                  << endl;*/
            this->WriteConfigurationInPDBFormat();
            this->SetTotalInteractionEnergy(this->CalculateTotalInteractionEnergy());
            this->SetlowestlRMSD(this->CalculatelowestlRMSD(v));
            v.push_back(*this);
            //cout << "\nID in configuration vector: " << v.size() -1 <<endl;
            found2 = true;
            found1 = true;
            found = true;
            try1 = 0;
            try2 = 0;
            return found;

        }//if they make a pair
        else{ //try to find second triangle
          try2--;

        }
      }//inner while
      try1--;
    }//outer while
    return found;
  }//end function


  //will be called from perturbation
  //c1 is to be perturbed
  //a copy of c1's parent is in dockedConfig (this) passed via this
  //dock within currnet triangle neighborhood of "this" and "c2"
  //dock "this" with "c2"
  void Configuration::DockAroundNeighborhood(std::vector <Triangle> &t1,
                                             std::vector <Triangle> &t2,
                                             int neighborTriangleId1,
                                             int neighborTriangleId2,
                                             Configuration& c1,
                                            Configuration& c2,
                                            std::vector <Configuration> &v){

    /*cout << "Inside DockAroundNeighborhood: ";
    cout << "before docking: ";
    cout << "curr reference and move unit triangle size: "<<endl;
    cout << this->GetTriangleInfoVectorSize() << endl;
    cout << c2.GetTriangleInfoVectorSize() << endl;*/
    //cout << "move and reference triangle id: ";
    //cout << c.GetReferenceTriangleId() << "\t" << c.GetMoveTriangleId()<<endl;

    real referenceRotationMatrix[9], moveRotationMatrix[9];
    Vector3 referenceTranslationVector3, moveTranslationVector3;

    //relative rotation and translation of reference frame wrt global
    t1[neighborTriangleId1].TransformTriangle
                (referenceRotationMatrix, referenceTranslationVector3);
    //relative rotation and translation of move frame wrt global
    t2[neighborTriangleId2].TransformTriangle
                (moveRotationMatrix, moveTranslationVector3);

    //calculate the final translation and rotation
    //final rotation
    //if two traingles are P(reference) and Q(move) wrt global
    //final rot from Q to P . R_PQ = R_WP.transpose() * R_WQ
    real referenceRotationMatrixTranspose[9];
    real finalRotationMatrix[9];
    //1. R_WP.transpose()
    TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
    //2. R_WP.transpose() * R_WQ
    MultMatrix(finalRotationMatrix,
               referenceRotationMatrixTranspose,
               moveRotationMatrix);

    //final translation
    //T = R_WP * A - R_WQ*D
    //T_PQ = R_WP.transpose * T
    Vector3 distanceOrigReferenceMove =
                    moveTranslationVector3 - referenceTranslationVector3;
    Vector3 finalTranslationVector3 =
    MultMatrixVector(referenceRotationMatrixTranspose,
                       distanceOrigReferenceMove);

    //transform move configuration on top of refernce configuration
    //using reference and move frame rotation
    //applying on move frame
    //return the whole docked configuration in dockedConfig object

    //initialize the docked config object with reference atom
    //information
    //replace this with copy vector

    unsigned long l = this->CopyTransformConfiguration
                (c2, finalRotationMatrix,
                 finalTranslationVector3);

    this->UpdateAfterDocking
        (c1.GetPdbId(), c1.GetAllChains(),
        c1.GetMoveChainId(), c1.GetReferenceChainId(),
        neighborTriangleId2, neighborTriangleId1,
        finalTranslationVector3, finalRotationMatrix,
        0, c1.referenceAtomEndIndex,
        c1.moveAtomStartIndex, c1.moveAtomEndIndex,
        c1.GetReferenceConfigurationInteractionEnergy());
        /*this->CalculateTotalInteractionEnergy(false,
                 this->GetTotalInteractionEnergy())*/
    /*cout << "Inside DockAroundNeighborhood: "
          << neighborTriangleId1 << "\t" << neighborTriangleId2 << endl;

    cout << "Updated value after docking: ";
    cout << this->GetInvolvedChainId() << "\t"
          << this->GetReferenceChainId() << "\t"
          << this->GetMoveChainId() << "\t"
          << this->GetReferenceTriangleId() << "\t"
          << this->GetMoveTriangleId() << "\t"
          << endl;
    cout << "size of the final structure: "
          << this->GetCurrentConfigurationSize() << endl;*/
    this->WriteConfigurationInPDBFormat();
    this->SetTotalInteractionEnergy(this->CalculateTotalInteractionEnergy());
    this->SetlowestlRMSD(this->CalculatelowestlRMSD(v));
    v.push_back(*this);
    //cout << "\nID in configuration vector: " << v.size() -1 <<endl;

  }


  int Configuration::GetDTScore(){
    return DTscore;
  }

  void Configuration::SetDTScore(int val_){
    DTscore = val_;
  }

	//idDock+ classifier
	int Configuration::IfTrueInteractionNew(){

		std::string str;
		real hydrophobicComp, acidicComp, hydrophilicComp, basicComp, conservationScore;
		real interfaceArea, interfaceAreaRatio;
		ofstream oFile;
		std::string fileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";

		//obtain interface properties
		std::string s1 = "python interfaceProperties.py";
		string s2 = s1 + " " + fileName + " " + GetInvolvedChainId();
		const char *cmd = s2.c_str();
		system(cmd);

		ifstream inFile;
		inFile.open("interfaceValues.dat");
		//first line is the properties name
		getline(inFile, str);
		//read the properties value
		inFile >> interfaceArea >> interfaceAreaRatio;
		inFile >> hydrophobicComp >> hydrophilicComp >> basicComp >> acidicComp;
		inFile >> conservationScore;
		//determine the true/false interaction
		std::vector<int> bc;
		bc.push_back(baseClassifier1(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier2(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier3(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier4(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier5(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier6(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier7(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier8(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier9(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));
		bc.push_back(baseClassifier10(interfaceArea, interfaceAreaRatio, hydrophobicComp, hydrophilicComp, basicComp, acidicComp, conservationScore));

		return CountOccurence(bc);

	}

  //check if true interaction
  //idDock classifier
  int Configuration::IfTrueInteractionOld(){

    std::string str;
    real hydrophobicComp, acidicComp, hydrophilicComp, basicComp, conservationScore;
    real interfaceArea, interfaceAreaRatio;
    ofstream oFile;
    std::string fileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb";

    //obtain interface properties
    std::string s1 = "python interfaceProperties.py";
    string s2 = s1 + " " + fileName + " " + GetInvolvedChainId();
    const char *cmd = s2.c_str();
    system(cmd);

    ifstream inFile;
    inFile.open("interfaceValues.dat");
    //first line is the properties name
    getline(inFile, str);
    //read the properties value
    inFile >> interfaceArea >> interfaceAreaRatio;
    inFile >> hydrophobicComp >> hydrophilicComp >> basicComp >> acidicComp;
    inFile >> conservationScore;
    //determine the true/false interaction
    if(hydrophobicComp > HYDROPHOBICTH){
     return 1;
    }
    if(hydrophobicComp <= HYDROPHOBICTH){
        if(acidicComp <= ACIDICTH){
            return 2;
        }
        else{
            if(conservationScore > JETTH){
                return 3;
            }
            else{
                if(interfaceAreaRatio > INTERFACERATIOTH){
                    return 0;
                }
                else{
                   if(hydrophilicComp <= HYDROPHILLICTH){
                     return 4;
                   }
                   else{
                       if(interfaceArea > INTERFACETH){
                            return 5;
                        }
                        else{
                            return 0;
                        }
                   }
                }

            }

        }

    }

}


  //dock 'this' with 'c' and get dockedConfig
  //only called from pairwise
  bool Configuration::Dock(Configuration& c,
                           Configuration& dockedConfig,
                           std::vector<Configuration> &v){

    //randomly sample two triangles from the triangle vector
    //cout << "Inside Dock: ";
    //cout << "before docking: ";
    //cout << "curr reference and move unit triangle size: "<<endl;
    //cout << this->GetTriangleInfoVectorSize() << endl;
    //cout << c.GetTriangleInfoVectorSize() << endl;

    int referenceTriangleId = IntUrandom(this->GetTriangleInfoVectorSize());
    int moveTriangleId = IntUrandom(c.GetTriangleInfoVectorSize());

    //find two geometrically fit triangle
    //if fit then transform
    int try1 = TRIANGLETRYCOUNT;
    int try2 = TRIANGLETRYCOUNT;

    bool found1 = false;
    bool found2 = false;
    bool found = false;

    while(try1 > 0  and found1 != true){

      int neighborTriangleId1 = findNeighborTriangle(referenceTriangleId, this->triangleInfoVector);
      if(neighborTriangleId1 < 0){
        try1--;
        //cout << "Not Found neighborid1 inside Dock\t" << try1 << endl;
        continue;

      }
      while(try2 > 0 and  found2 != true){
        int neighborTriangleId2 = findNeighborTriangle(moveTriangleId, c.triangleInfoVector);
        if(neighborTriangleId2 < 0){
          try2--;
          //cout << "Not Found neighborid2 inside Dock\t" << try2 << endl;
          continue;

        }
        //if geometric shape complement and other geometric constraints meet
        if(triangleInfoVector[neighborTriangleId1].CheckShapeComplementarity
          (c.triangleInfoVector[neighborTriangleId2])
          &&  (triangleInfoVector[neighborTriangleId1].PairwiseTriangle
            (c.triangleInfoVector[neighborTriangleId2]))){
            //if shape and geometrically complementary
            //transform move triangle on top of reference
            //obtain final rotation and translation
            real referenceRotationMatrix[9], moveRotationMatrix[9];
            Vector3 referenceTranslationVector3, moveTranslationVector3;

            //relative rotation and translation of reference frame wrt global
            this->triangleInfoVector[neighborTriangleId1].TransformTriangle
                  (referenceRotationMatrix, referenceTranslationVector3);
            //relative rotation and translation of move frame wrt global
            c.triangleInfoVector[neighborTriangleId2].TransformTriangle
                  (moveRotationMatrix, moveTranslationVector3);
            //update triangle info
            this->SetReferenceTriangleId(neighborTriangleId1);
            c.SetMoveTriangleId(neighborTriangleId2);

            //calculate the final translation and rotation
            //final rotation
            //if two traingles are P(reference) and Q(move) wrt global
            //final rot from Q to P . R_PQ = R_WP.transpose() * R_WQ
            real referenceRotationMatrixTranspose[9];
            real finalRotationMatrix[9];
            //1. R_WP.transpose()
            TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
            //2. R_WP.transpose() * R_WQ
            MultMatrix(finalRotationMatrix,
                       referenceRotationMatrixTranspose,
                       moveRotationMatrix);

            //final translation
            //T = R_WP * A - R_WQ*D
            //T_PQ = R_WP.transpose * T
            Vector3 distanceOrigReferenceMove =
                          moveTranslationVector3 - referenceTranslationVector3;
            Vector3 finalTranslationVector3 =
            MultMatrixVector(referenceRotationMatrixTranspose,
                             distanceOrigReferenceMove);

            //transform move configuration on top of refernce configuration
            //using reference and move frame rotation
            //applying on move frame
            //return the whole docked configuration in dockedConfig object

            //initialize the docked config object with reference and move atom
            //information
            //replace this with copy vector
            for(int k = 0; k < this->GetCurrentConfigurationSize(); k++){
                Atom tempAtom
                   (this->currentConfigurationInfo[k].GetAtomSerialNumber(),
                    this->currentConfigurationInfo[k].GetAtomName(),
                    this->currentConfigurationInfo[k].GetAltLoc(),
                    this->currentConfigurationInfo[k].GetResName(),
                    this->currentConfigurationInfo[k].GetChainID(),
                    this->currentConfigurationInfo[k].GetResSeq(),
                    this->currentConfigurationInfo[k].GetX(),
                    this->currentConfigurationInfo[k].GetY(),
                    this->currentConfigurationInfo[k].GetZ(),
                    this->currentConfigurationInfo[k].GetOccupancy(),
                    this->currentConfigurationInfo[k].GetTempFactor(),
                    this->currentConfigurationInfo[k].GetElement(),
                    false);
                dockedConfig.currentConfigurationInfo.push_back(tempAtom);
            }

            unsigned long l = this->TransformConfiguration
                  (c, finalRotationMatrix,
                   finalTranslationVector3, dockedConfig);

            dockedConfig.UpdateAfterDocking
                  (this->GetPdbId(), dockedConfig.GetAllChains(),
                  c.GetInvolvedChainId(), this->GetInvolvedChainId(),
                  neighborTriangleId2, neighborTriangleId1,
                  finalTranslationVector3, finalRotationMatrix,
                  0, this->GetCurrentConfigurationSize()-1,
                   this->GetCurrentConfigurationSize(),
                   dockedConfig.GetCurrentConfigurationSize()-1,
                   this->GetReferenceConfigurationInteractionEnergy());
                   /*dockedConfig.CalculateTotalInteractionEnergy(true,
                    this->GetTotalInteractionEnergy())*/
            /*cout << "NeighborTriangleIds: "
                  << neighborTriangleId1 << "\t" << neighborTriangleId2 << endl;
            cout << "Updated value after docking: ";
            cout << dockedConfig.GetInvolvedChainId() << "\t"
                  //<< dockedConfig.GetReferenceTriangleId() << "\t"
                  << dockedConfig.GetCurrentConfigurationSize() << "\t"
                  << endl;*/
            dockedConfig.WriteConfigurationInPDBFormat();
            if(DT==1){
              int ifTrue = dockedConfig.IfTrueInteractionNew();
              if(ifTrue){

                if(CLUSTER){ //if cluster = on then calc e for only leaders later
                  dockedConfig.SetDTScore(ifTrue);
                  dockedConfig.SetTotalInteractionEnergy(0.0);
                  dockedConfig.SetlowestlRMSD(0.0);

                }
                else{ // else do as usual no clustering
                  //dockedConfig.WriteConfigurationInPDBFormat();
                  if(PAIRWISEENERGY==0){//only DT no energy
                    dockedConfig.SetDTScore(ifTrue);
                    dockedConfig.SetTotalInteractionEnergy(0.0); //not measuring energy
                  }
                  else{ //energy on just regular options calc e for all seeds
                    dockedConfig.SetDTScore(ifTrue);
                    dockedConfig.SetTotalInteractionEnergy(
                              dockedConfig.CalculateTotalInteractionEnergy());
					dockedConfig.SetlowestlRMSD(dockedConfig.CalculatelowestlRMSD(v));
                  }
                  //dockedConfig.SetlowestlRMSD(dockedConfig.CalculatelowestlRMSD(v));
                }
                v.push_back(dockedConfig);
                //cout << "\ntrue interaction found ID in configuration vector: "
                //<< v.size() -1 <<endl;
                //cout << "\nDT score: "
                //<< dockedConfig.GetDTScore() <<endl;
                found2 = true;
                found1 = true;
                found = true;
                try1 = 0;
                try2 = 0;
                return found;
              }
              else{//not a true interaction //sample another triangle2
                try2--;
                dockedConfig.currentConfigurationInfo.clear();
                //cout << "Not a true interaction\n";

              }//else of ifTrue
            }//DT on
            else{ //DT off
              dockedConfig.SetTotalInteractionEnergy(dockedConfig.CalculateTotalInteractionEnergy());
              dockedConfig.SetlowestlRMSD(dockedConfig.CalculatelowestlRMSD(v));
              v.push_back(dockedConfig);
              //cout << "\nDT off: so accept everything\n";
              //cout << "ID in configuration vector: " << v.size() -1 <<endl;
              found2 = true;
              found1 = true;
              found = true;
              try1 = 0;
              try2 = 0;
              return found;
            }//DT off

        }//if they make a pair
        else{ //try to find second triangle
          try2--;

        }
      }//inner while
      try1--;
    }//outer while
    return found;
  }//end function


  //misc utilty functions
  //copy src from start to end to this
  void Configuration::CopyConfigurationInfo(std::vector <Atom> &src,
                                            unsigned long start,
                                            unsigned long end){
      for (int i = start; i <= end; i++){
        this->currentConfigurationInfo.push_back(src[i]);
      }

  }

  /*void Configuration::SetReferenceConfigurationIndex
  (unsigned long refConfigIndex_){
    referenceConfigurationIndex = refConfigIndex_;

  }

  unsigned long GetReferenceConfigurationIndex(){
    return referenceConfigurationIndex;
  }*/

  //for d2d if genNo == 1000 then
  //if j >=14397 and j <= 21954 print C
  //if j>21954 and j <= 28793 print D
  void Configuration::WriteConfigurationInPDBFormatInLogFile
                            (std::string fileName,
                             int generationNo, int indvNo){

    ofstream oFile;
    oFile.open(StringToChar(fileName));
    string spacer(" ");
    int dash;
    oFile << "Generation " << generationNo << " Individual " << indvNo
          << " Energy " << GetTotalInteractionEnergy()
          << " RMSD " << GetlowestlRMSD() << endl;

    for(int j = 0 ; j < GetCurrentConfigurationSize();j++){
      stringstream atom_stream;
      if (currentConfigurationInfo[j].GetAtomName().size() > 3){
          dash = 0;
      }
      else{
        dash = 1;
      }
      atom_stream << setiosflags(ios::left) << (spacer.c_str()+1)-dash
                  << setw(4-dash) << currentConfigurationInfo[j].GetAtomName()
                  << " ";
      oFile << setiosflags(ios:: right)
            << setiosflags(ios::fixed) << setprecision(THREE);
      oFile << "ATOM  " << setw(FIVE)
            << currentConfigurationInfo[j].GetAtomSerialNumber() <<  " "
            << atom_stream.str() << currentConfigurationInfo[j].GetResName();

      oFile << " " << currentConfigurationInfo[j].GetChainID() << " ";
      oFile  << setw(THREE)<< currentConfigurationInfo[j].GetResSeq() << "    ";
      oFile << setw(EIGHT) << currentConfigurationInfo[j].GetX()
            << setw(EIGHT) << currentConfigurationInfo[j].GetY()
            << setw(EIGHT) << currentConfigurationInfo[j].GetZ();

      oFile << "  " << "1.00" << endl;
    }
    oFile << "TER" << endl;
    //done generating a temp.pdb file
    oFile.close();
  }







//  bool Configuration::RandomDock(Configuration *dockedConfig,
//                                 std::vector <Configuration*>&v){
//
//    //convert finalrotation to elt/quat
//    //quat to randomquat
//    //randomquat to randommatrix
//    //transpose finalrotation and multiply with randommatrix and finalrandomrotmatrix
//    //final translation to randomtranslation
//    //diff between finaltranslation - randomtranslation
//    //multiply with randomrotation with diff and finalrandomtranslation
//    //cout << "Inside RandomDock: " << endl;
//    real resElt[4], randomizeQuat[4];
//    std::vector<real> resQuat;
//    real randomRotationMatrix[9];
//    real rotationMatrixTranspose[9];
//    real finalRandomRotationMatrix[9];
//    real randomTranslation[3];
//    //unsigned long moveIndex = this->GetMoveTriangleId();
//    //obtain the point
//
//    //rotation matrix to quat
//    ConvertFromMatrixToQuat(resQuat, rotationMatrix);
//    //quat to elt
//    resElt[0] = resQuat[1];
//    resElt[1] = resQuat[2];
//    resElt[2] = resQuat[3];
//    resElt[3] = resQuat[0];
//
//    //randomize quat
//    QuaternionRandomNeigh( resElt,  DTHETA, DPHI, randomizeQuat );
//    ConvertFromEltToMatrix(randomRotationMatrix, randomizeQuat);
//
//    //TransposeMatrix(rotationMatrixTranspose, rotationMatrix);
//    //MultMatrix(finalRandomRotationMatrix,
//               //rotationMatrixTranspose, randomRotationMatrix);
//
//    PositionRandomNeigh(translationVector, TRANSLATIONTHRESHOLD,
//                      TRANSLATIONTHRESHOLD, TRANSLATIONTHRESHOLD,
//                      randomTranslation);
//
//
//    Vector3 randomTranslationVector3(randomTranslation[0],
//                                     randomTranslation[1],
//                                     randomTranslation[2]);
//
//
//    Vector3 diffTranslation = randomTranslationVector3 - translationVector;
//    //Vector3 finalRandomTranslationVector3 =
//                //MultMatrixVector(finalRandomRotationMatrix, diffTranslation);
//
//
//    //Vector3 finalRandomTranslationVector3 = randomTranslationVector3;
//    dockedConfig.CopyConfigurationInfo(this->currentConfigurationInfo,
//                                        this->referenceAtomStartIndex,
//                                        this->referenceAtomEndIndex);
//
//    //for(int k = 0 ;  k < 9; k++)randomRotationMatrix[k] = 1.0;
//    unsigned long l = this->TransformConfigurationWithIndex
//                  (randomRotationMatrix,
//                   diffTranslation, dockedConfig,
//                   this->moveAtomStartIndex, this->moveAtomEndIndex);
//    //cout << "after transformation size of transformed: " << l << endl;
//
//    dockedConfig.referenceAtomStartIndex = this->referenceAtomStartIndex;
//    dockedConfig.referenceAtomEndIndex = this->referenceAtomEndIndex;
//    dockedConfig.moveAtomStartIndex = this->moveAtomStartIndex;
//    dockedConfig.moveAtomEndIndex = this->moveAtomEndIndex;
//
//    /*cout <<"involved chains: "<< dockedConfig.GetInvolvedChainId()<< endl;
//    cout << dockedConfig.referenceAtomStartIndex << endl;
//    cout << dockedConfig.referenceAtomEndIndex << endl;
//    cout << dockedConfig.moveAtomStartIndex << endl;
//    cout << dockedConfig.moveAtomEndIndex << endl;
//    */
//    dockedConfig.UpdateAfterDocking
//          (this->GetPdbId(), this->GetAllChains(),
//          this->GetMoveChainId(), this->GetReferenceChainId(),
//          this->GetMoveTriangleId(), this->GetReferenceTriangleId(),
//          randomTranslationVector3, randomRotationMatrix,
//          this->referenceAtomStartIndex, this->referenceAtomEndIndex,
//          this->moveAtomStartIndex, this->moveAtomEndIndex,
//          this->GetReferenceConfigurationInteractionEnergy()
//          /*this->CalculateTotalInteractionEnergy(false,
//                   this->GetTotalInteractionEnergy())*/);
//
//    //cout << "after docking: chains" << dockedConfig.GetInvolvedChainId() << endl;
//    //cout << "after docking: size" << dockedConfig.GetCurrentConfigurationSize() << endl;
//
//    dockedConfig.WriteConfigurationInPDBFormat();
//    dockedConfig.SetTotalInteractionEnergy(dockedConfig.CalculateTotalInteractionEnergy());
//    dockedConfig.SetlowestlRMSD(dockedConfig.CalculatelowestlRMSD(v));
//    v.push_back(dockedConfig);
//
//    /*cout << "energy in modify: " << dockedConfig.GetTotalInteractionEnergy();
//    cout << "rmsd after docking: " << dockedConfig.CalculatelowestlRMSD(v);
//	*/
//
//   }


  /*unsigned long Configuration::GetParentId(){
    return parentId;
  }
  void Configuration::SetParentId(unsigned long parentId_){
    parentId = parentId_;

  }*/

/*bool Configuration::RandomDock(Configuration *dockedConfig,
                                 std::vector <Configuration*>&v){

    //calculate center of mass
    //move to origin
    //rotate
    //translate back
    std::vector<Atom*> a;
    for(int i = moveAtomStartIndex ; i <=moveAtomEndIndex; i++ ){
      a.push_back(this->currentConfigurationInfo[i]);
    }
    Vector3 com = CalculateCenterOfMass(a);





    dockedConfig.CopyConfigurationInfo(this->currentConfigurationInfo,
                                        this->referenceAtomStartIndex,
                                        this->referenceAtomEndIndex);

    //for(int k = 0 ;  k < 9; k++)randomRotationMatrix[k] = 1.0;
    unsigned long l = this->TransformConfigurationWithIndex
                  (randomRotationMatrix,
                   diffTranslation, dockedConfig,
                   this->moveAtomStartIndex, this->moveAtomEndIndex);


   }*/

   /* ******************** codes for randomDock *********************/


    //COM
    Vector3 Configuration::CalculateCenterOfMass(std::vector<Atom>& a){

      real centerOfMassX = 0.0;
      real centerOfMassY = 0.0;
      real centerOfMassZ = 0.0;
      unsigned long l = a.size();
      for(int i =0 ; i < l; i++){
        centerOfMassX += a[i].GetX();
        centerOfMassY += a[i].GetY();
        centerOfMassZ += a[i].GetZ();

      }
      Vector3 centerOfMass(centerOfMassX/l, centerOfMassY/l, centerOfMassZ/l);
      return centerOfMass;

    }
    //c2: moving unit
    //docked Config: computed configuration
    //v: configuration vector
    void Configuration::RandomDock(Configuration& c2,
                                   Configuration& dockedConfig,
                                   std::vector <Configuration>& v){

      //cout << "Inside RandomDock\n";
      //push the this configuration in to dockedConfig
      //push the reference config in dockedConfig
      std::vector <Atom> tempAtomVect;
      this->GetCurrentConfigurationInfo(tempAtomVect);
      dockedConfig.SetCurrentConfigurationInfo(tempAtomVect);
      //cout << tempAtomVect.size() << endl;
      //cout << dockedConfig.GetCurrentConfigurationSize() << endl;

      //randomly rotate
      //move it back
      //random translate
      //copy move unit configuration in Atom vector a
      std::vector<Atom> a;
      for(int i = 0;
          i < c2.GetCurrentConfigurationSize();
          i++ ){
          a.push_back(c2.currentConfigurationInfo[i]);
      }
      //cout << a.size() << endl;
      //calculate center of mass of vector a
      Vector3 com = CalculateCenterOfMass(a);
      //cout << "COM " << com << endl;
      //move a to origin
      Vector3 origin(-1.0 * com[0], -1.0 * com[1], -1.0 * com[2]);
      //cout << "origin " << origin << endl;
      //cout << "Move to Origin\n";
      for(int i = 0 ; i < a.size(); i++){
        a[i].SetX(a[i].GetX()+origin[0]);
        a[i].SetY(a[i].GetY()+origin[1]);
        a[i].SetZ(a[i].GetZ()+origin[2]);
        //cout << a[i].GetX() << "\t";
        //cout << a[i].GetY() << "\t";
        //cout << a[i].GetZ() << "\t";
      }
      //cout << endl;

      //rotate slightly along X, Y and Z axis
      real degree = RealRandRange(-1.0 * RANDDEG, RANDDEG);
      //cout << "degree " << degree << endl;
      for(int i = 0;i < a.size(); i++){
        //cout << a[i].GetX() << "\t";
        //cout << a[i].GetY() << "\t";
        //cout << a[i].GetZ() << "\t";
        real pz[3];
        a[i].RotateAlongAxis(pz, degree, a[i].GetX(), a[i].GetY(),
                                            a[i].GetZ());
        //cout << pz[0] << "\t";
        //cout << pz[1] << "\t";
        //cout << pz[2] << "\t";
        a[i].SetX(pz[0]);
        a[i].SetY(pz[1]);
        a[i].SetZ(pz[2]);
        //cout << a[i].GetX() << "\t";
        //cout << a[i].GetY() << "\t";
        //cout << a[i].GetZ() << "\t";
      }
      //cout << endl;

      //move back again from origin to previous position
      //cout << "move back from origin\n";
      for(int i = 0 ; i < a.size(); i++){
        a[i].SetX(a[i].GetX()+com[0]);
        a[i].SetY(a[i].GetY()+com[1]);
        a[i].SetZ(a[i].GetZ()+com[2]);
        //cout << a[i] << "\t";
      }
      //cout << endl;
      //translate later
      //push back this to the docked Config
      for (int i = 0 ;i < a.size(); i++){
        dockedConfig.currentConfigurationInfo.push_back(a[i]);
      }
      //rotation and traslation is not needed anywhere
      //putting a place holded
      //need to change
      real finalRotationMatrix[9] = {0,0,0,
                                    0,0,0,
                                    0, 0, 0};
      Vector3 finalTranslationVector(0,0,0);
      //write a different update function
      dockedConfig.UpdateAfterDocking(dockedConfig.GetPdbId(),
                                       dockedConfig.GetAllChains(),
                                       dockedConfig.GetMoveChainId(),
                                       dockedConfig.GetReferenceChainId(),
                                       0, 0,
                                       finalTranslationVector, finalRotationMatrix,
                                       0,
                                       this->GetCurrentConfigurationSize()-1,
                                       this->GetCurrentConfigurationSize(),
                                       dockedConfig.GetCurrentConfigurationSize()-1,
                                       0.0);

      //dockedConfig.SetReferenceChainId(this->GetInvolvedChainId());
      //dockedConfig.SetMoveChainId(c2.GetInvolvedChainId());

      //write pdb for foldx energy calculateion
      dockedConfig.WriteConfigurationInPDBFormat();
      //calaute foldX energy
      dockedConfig.SetTotalInteractionEnergy(dockedConfig.CalculateTotalInteractionEnergy());
      //calucluate lRMSD form native
      dockedConfig.SetlowestlRMSD(dockedConfig.CalculatelowestlRMSD(v));
      v.push_back(dockedConfig);
//      cout <<"Updated Value after docking:\n";
//      cout << dockedConfig.GetInvolvedChainId() << endl;
//      cout << dockedConfig.GetCurrentConfigurationSize() << endl;
//      cout << dockedConfig.GetlowestlRMSD() << endl;
//      cout << dockedConfig.GetTotalInteractionEnergy() << endl;
//      cout << v.size() << endl;
   }

    void Configuration::PrintAConfiguration(){
      std::vector <Atom>a;
      cout << "\n\nPrinting Configurations:\n";
      this->GetCurrentConfigurationInfo(a);
      for(int i = 0 ; i < a.size(); i++){
        cout << a[i].GetX() << "\t";
        cout << a[i].GetY() << "\t";
        cout << a[i].GetZ() << "\t";

      }
      cout << "\n END PRINTING ATOMS\n";

    }//end of function

    //clustering variables
    int Configuration::GetMembershipStatus(){
      return membershipStatus;

    } //0 = not selected, 1 = leader, 2 = memeber
    void Configuration::SetMembershipStatus(int status_){

      membershipStatus = status_;

    } //0 = not selected, 1 = leader, 2 = memeber

    const unsigned long Configuration::GetClusterSize() const{

      return clusterSize;

    } //applicable for leaders only

    void Configuration::SetClusterSize(int size_){
      clusterSize = size_;

    }

    void Configuration::IncreaseClusterSize(void){
      clusterSize++;

    }


    /* ********************* d2d stuff ************************ */
    real Configuration::GetTotalTimingInfo(){
      return timingInfo;

    }

    void Configuration::SetTotalTimingInfo(real time_){
      timingInfo = time_;

    }

    real Configuration::GetEnergyTimingInfo(){
      return energyTimingInfo;

    }
    void Configuration::SetEnergyTimingInfo(real time_){
      energyTimingInfo = time_;

    }



    void Configuration::PrepareForDockingInSpecificRegion(){
      WriteConfigurationInPDBFormat();
      CalculateConnollyRepresentation();
      CalculateSpecificShuoRepresentation();
      CalculateTriangleRepresentation();

    }
	//sampling is only done from the specific TM regions
    void Configuration::CalculateSpecificShuoRepresentation(){
      int TM1[] = {596, 946};
      int TM2[] = {1110, 1483};
      int TM3[] = {1731, 2067};
      int TM4[] = {2428, 2772};
      int TM5[] = {2991, 3427};
      int TM6[] = {6044, 6420};
      int TM7[] = {6631, 6980};

      //calculate shuo representation
      CalculateShuoRepresentation();
      //I will modify the file with TM points and then
      //rename as the same name as shuo file
      //then pass it to the triangle generation stage
      //filename Format:  shuoFileName = GetPdbId() + "_" + GetInvolvedChainId() + ".pdb.shuo";
      std::string surfaceFileName = GetPdbId() + "_"
                              + GetInvolvedChainId()
                              + ".pdb.shuo";
      std::vector <Surface> surfaceVector, temp;
      ReadShuo(surfaceFileName, surfaceVector);
      RemoveFile(surfaceFileName);
      ofstream sFile;
      sFile.open(StringToChar(surfaceFileName));

      //process surface vector
      for(int i = 0; i < surfaceVector.size(); i++){
        if(
        ((surfaceVector[i].atom1 >= TM1[0] and surfaceVector[i].atom1 <= TM1[1])
      or (surfaceVector[i].atom2 >= TM1[0] and surfaceVector[i].atom2 <= TM1[1])
      or (surfaceVector[i].atom3 >= TM1[0] and surfaceVector[i].atom3 <= TM1[1]))
      or

      ((surfaceVector[i].atom1 >= TM2[0] and surfaceVector[i].atom1 <= TM2[1])
      or (surfaceVector[i].atom2 >= TM2[0] and surfaceVector[i].atom2 <= TM2[1])
      or (surfaceVector[i].atom3 >= TM2[0] and surfaceVector[i].atom3 <= TM2[1]))
      or

      ((surfaceVector[i].atom1 >= TM3[0] and surfaceVector[i].atom1 <= TM3[1])
      or (surfaceVector[i].atom2 >= TM3[0] and surfaceVector[i].atom2 <= TM3[1])
      or (surfaceVector[i].atom3 >= TM3[0] and surfaceVector[i].atom3 <= TM3[1]))

      or

        ((surfaceVector[i].atom1 >= TM4[0] and surfaceVector[i].atom1 <= TM4[1])
      or (surfaceVector[i].atom2 >= TM4[0] and surfaceVector[i].atom2 <= TM4[1])
      or (surfaceVector[i].atom3 >= TM4[0] and surfaceVector[i].atom3 <= TM4[1]))

      or

      ((surfaceVector[i].atom1 >= TM5[0] and surfaceVector[i].atom1 <= TM5[1])
      or (surfaceVector[i].atom2 >= TM5[0] and surfaceVector[i].atom2 <= TM5[1])
      or (surfaceVector[i].atom3 >= TM5[0] and surfaceVector[i].atom3 <= TM5[1]))

      or

      ((surfaceVector[i].atom1 >= TM6[0] and surfaceVector[i].atom1 <= TM6[1])
      or (surfaceVector[i].atom2 >= TM6[0] and surfaceVector[i].atom2 <= TM6[1])
      or (surfaceVector[i].atom3 >= TM6[0] and surfaceVector[i].atom3 <= TM6[1]))

      or
      ((surfaceVector[i].atom1 >= TM7[0] and surfaceVector[i].atom1 <= TM7[1])
      or (surfaceVector[i].atom2 >= TM7[0] and surfaceVector[i].atom2 <= TM7[1])
      or (surfaceVector[i].atom3 >= TM7[0] and surfaceVector[i].atom3 <= TM7[1]))
      ){
          sFile << surfaceVector[i].atom1 << "\t";
          sFile << surfaceVector[i].atom2 << "\t";
          sFile << surfaceVector[i].atom3 << "\t";
          sFile << surfaceVector[i].x << "\t";
          sFile << surfaceVector[i].y << "\t";
          sFile << surfaceVector[i].z << "\n";


        }//if

      }//for
      //write modified shuo information in shuoFileName

      sFile.close();
    }

    //dock 'this' with 'c' and get dockedConfig
    bool Configuration::d2dDock(Configuration& c,
                    Configuration& dockedConfig,
                    std::vector <Configuration> &v, int totalConfigs){

      //randomly sample two triangles from the triangle vector
      cout << "Inside d2d Dock: ";
      MyTimer runtime;
      runtime.start();
      int referenceTriangleId = IntUrandom(this->GetTriangleInfoVectorSize());
      int moveTriangleId = IntUrandom(c.GetTriangleInfoVectorSize());

      //find two geometrically fit triangle
      //if fit then transform
      int try1 = TRIANGLETRYCOUNT;
      int try2 = TRIANGLETRYCOUNT;

      bool found1 = false;
      bool found2 = false;
      bool found = false;

      while(try1 > 0  and found1 != true){

        int neighborTriangleId1 = findNeighborTriangle(referenceTriangleId, this->triangleInfoVector);
        if(neighborTriangleId1 < 0){
          try1--;
          cout << "Not Found neighborid1 inside Dock\t" << try1 << endl;
          continue;

        }
        while(try2 > 0 and  found2 != true){
          int neighborTriangleId2 = findNeighborTriangle(moveTriangleId, c.triangleInfoVector);
          if(neighborTriangleId2 < 0){
            try2--;
            cout << "Not Found neighborid2 inside Dock\t" << try2 << endl;
            continue;

          }
          //if geometric shape complement and other geometric constraints meet
          if(triangleInfoVector[neighborTriangleId1].CheckShapeComplementarity
            (c.triangleInfoVector[neighborTriangleId2])
            &&  (triangleInfoVector[neighborTriangleId1].PairwiseTriangle
              (c.triangleInfoVector[neighborTriangleId2]))){
              //if shape and geometrically complementary
              //transform move triangle on top of reference
              //obtain final rotation and translation
              real referenceRotationMatrix[9], moveRotationMatrix[9];
              Vector3 referenceTranslationVector3, moveTranslationVector3;

              //relative rotation and translation of reference frame wrt global
              this->triangleInfoVector[neighborTriangleId1].TransformTriangle
                    (referenceRotationMatrix, referenceTranslationVector3);
              //relative rotation and translation of move frame wrt global
              c.triangleInfoVector[neighborTriangleId2].TransformTriangle
                    (moveRotationMatrix, moveTranslationVector3);
              //update triangle info
              this->SetReferenceTriangleId(neighborTriangleId1);
              c.SetMoveTriangleId(neighborTriangleId2);

              //calculate the final translation and rotation
              //final rotation
              //if two traingles are P(reference) and Q(move) wrt global
              //final rot from Q to P . R_PQ = R_WP.transpose() * R_WQ
              real referenceRotationMatrixTranspose[9];
              real finalRotationMatrix[9];
              //1. R_WP.transpose()
              TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
              //2. R_WP.transpose() * R_WQ
              MultMatrix(finalRotationMatrix,
                         referenceRotationMatrixTranspose,
                         moveRotationMatrix);

              //final translation
              //T = R_WP * A - R_WQ*D
              //T_PQ = R_WP.transpose * T
              Vector3 distanceOrigReferenceMove =
                            moveTranslationVector3 - referenceTranslationVector3;
              Vector3 finalTranslationVector3 =
              MultMatrixVector(referenceRotationMatrixTranspose,
                               distanceOrigReferenceMove);

              //transform move configuration on top of refernce configuration
              //using reference and move frame rotation
              //applying on move frame
              //return the whole docked configuration in dockedConfig object

              //initialize the docked config object with reference and move atom
              //information
              //replace this with copy vector
              for(int k = 0; k < this->GetCurrentConfigurationSize(); k++){
                  Atom tempAtom
                     (this->currentConfigurationInfo[k].GetAtomSerialNumber(),
                      this->currentConfigurationInfo[k].GetAtomName(),
                      this->currentConfigurationInfo[k].GetAltLoc(),
                      this->currentConfigurationInfo[k].GetResName(),
                      this->currentConfigurationInfo[k].GetChainID(),
                      this->currentConfigurationInfo[k].GetResSeq(),
                      this->currentConfigurationInfo[k].GetX(),
                      this->currentConfigurationInfo[k].GetY(),
                      this->currentConfigurationInfo[k].GetZ(),
                      this->currentConfigurationInfo[k].GetOccupancy(),
                      this->currentConfigurationInfo[k].GetTempFactor(),
                      this->currentConfigurationInfo[k].GetElement(),
                      false);
                  dockedConfig.currentConfigurationInfo.push_back(tempAtom);
              }

              unsigned long l = this->TransformConfiguration
                    (c, finalRotationMatrix,
                     finalTranslationVector3, dockedConfig);

              dockedConfig.UpdateAfterDocking
                    (this->GetPdbId(), dockedConfig.GetAllChains(),
                    c.GetInvolvedChainId(), this->GetInvolvedChainId(),
                    neighborTriangleId2, neighborTriangleId1,
                    finalTranslationVector3, finalRotationMatrix,
                    0, this->GetCurrentConfigurationSize()-1,
                     this->GetCurrentConfigurationSize(),
                     dockedConfig.GetCurrentConfigurationSize()-1,
                     this->GetReferenceConfigurationInteractionEnergy());
                     /*dockedConfig.CalculateTotalInteractionEnergy(true,
                      this->GetTotalInteractionEnergy())*/
              cout << "NeighborTriangleIds: "
                    << neighborTriangleId1 << "\t" << neighborTriangleId2 << endl;
              cout << "Updated value after docking: ";
              cout << dockedConfig.GetInvolvedChainId() << "\t"
                    //<< dockedConfig.GetReferenceTriangleId() << "\t"
                    << dockedConfig.GetCurrentConfigurationSize() << "\t"
                    << endl;
              dockedConfig.WriteConfigurationInPDBFormat();
              clock_t energyTimeStart = clock_t();
              dockedConfig.SetTotalInteractionEnergy(
                                dockedConfig.CalculateTotalInteractionEnergy());

              clock_t energyTimeEnd = clock_t();
              dockedConfig.SetEnergyTimingInfo(
                (energyTimeEnd-energyTimeStart)/((double)CLOCKS_PER_SEC));

              v.push_back(dockedConfig);
              cout << "\nID in configuration vector: "
              << v.size() -1 <<endl;
              found2 = true;
              found1 = true;
              found = true;
              try1 = 0;
              try2 = 0;
              return found;

          }//if they make a pair
          else{ //try to find second triangle
            try2--;

          }
        }//inner while
        try1--;
      }//outer while
      return found;
  }//end function

}//end namespace

