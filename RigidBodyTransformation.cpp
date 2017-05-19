/*

Author: Irina Hashmi (ihashmi@masonlive.gmu.edu)
Advisor: Prof. Amarda Shehu
Time and place: Sep, 2014, Department of Computer Science, George Mason University, VA, USA
Goal:
Note:
Compile:
Run:
Copyright (c) 2014 Irina Hashmi. All rights reserved.

*/


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
    this->SetTotalInteractionEnergy(this->CalculateTotalInteractionEnergy());
    this->SetlowestlRMSD(this->CalculatelowestlRMSD(v));
    v.push_back(*this);


}

