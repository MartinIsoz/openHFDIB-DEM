/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    class for addition and removal of immersedBodies

SourceFiles
    addModelDistribution.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelDistribution_H
#define addModelDistribution_H

#include "addModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class addModelDistribution Declaration
\*---------------------------------------------------------------------------*/

class addModelDistribution
:
    public addModel
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar>  DynamicScalarList;
    const dictionary&           addModelDict_;
    word                        addMode_;
    bool                        bodyAdded_;
    autoPtr<geomModel> geomModel_;

    IOdictionary           distributionDict_;
    const scalarList       distribution_;
    const scalarList       particleSize_;
    const scalar       convertToMeters_;

    DynamicScalarList  addedParticlesSize_;
    scalar               volumeOfAddedBodies_;

    const dictionary&   coeffsDict_;
    const scalar       stlBaseSize_;

	word			    addDomain_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& addModeICoeffs_;

    label       useNTimes_;
    scalar      timeBetweenUsage_;
    label		partPerAdd_;
	scalar		fieldValue_;
    label		addedOnTimeLevel_;
    label		partPerAddTemp_;

	word		zoneName_;
    vector		minBound_;
    vector		maxBound_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		succesfulladition_;
	bool		restartPartCountTemp_;
	bool		reapeatedAddition_;
	bool		firstTimeRunning_;
	bool 		cellZoneActive_;
	bool		boundBoxActive_;
	Field<label>octreeField_;
	bool		timeBased_;
	bool		fieldBased_;
	scalar		fieldCurrentValue_;

    List<DynamicLabelList>  cellsInBoundBox_;
    List<pointField>        cellZonePoints_;
    boundBox                cellZoneBounds_;

    bool        allActiveCellsInMesh_;

    label       nGeometricD_;
    Vector<label>      geometricD_;

    Random          randGen_;

    protected:


    public:

    // Constructors ------------------------------------------------
    addModelDistribution
    (
        const dictionary& addModelDict,
        const Foam::dynamicFvMesh& mesh,
        const Vector<label> geomDir,
        geomModel* bodyGeomModel
    );

    // Destructors -------------------------------------------------
    ~addModelDistribution();

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    geomModel* addBody(const volScalarField& body);

    Tuple2<label, scalar> returnScaleFactor();

	vector computeCoM(const triSurfaceMesh& bodySurfMesh);

	bool isBodyInBoundBox
    (
        vector minBound_,
        vector maxBound_,
        const triSurfaceMesh& bodySurfMesh
    );

	void initializeCellZone();

	void updateCellZoneBoundBox();

    void initializeBoundBox();

    void recreateBoundBox();

	labelList getBBoxCellsByOctTree
    (
        label cellToCheck,
        bool& insideBB,
        vector& bBoxMin,
        vector& bBoxMax,
        List<DynamicLabelList>& bBoxCells
    );

	scalar checkLambdaFraction(const volScalarField& body);

    bool checkIfAdd() {return true;}

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

    scalar returnRandomAngle();

    vector returnRandomRotationAxis();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
