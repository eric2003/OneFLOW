#include "EulerDomain.h"
#include "AccelViews.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace
{

using namespace ONEFLOW;

EulerDomainProblem Problem()
{
    EulerDomainProblem result;
    result.nCells = 32;
    result.dt = 0.0005;
    result.dx = 1.0 / result.nCells;
    return result;
}

EulerDomainProblem FiveEquationProblem()
{
    EulerDomainProblem result = Problem();
    result.nEquations = 5;
    return result;
}

TEST( EulerDomainContract, ValidatesBackendNeutralProblemAndViews )
{
    const EulerDomainProblem problem = Problem();
    Real values[ 3 * 32 ] = {};
    EulerDomainConstFieldView input{ 32, 3, values };
    EulerDomainFieldView output{ 32, 3, values };

    EXPECT_NO_THROW( ValidateEulerDomainProblem( problem ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, input ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, output ) );
}

TEST( EulerDomainContract, ValidatesFiveEquationProblem )
{
    const EulerDomainProblem problem = FiveEquationProblem();
    Real values[ 5 * 32 ] = {};
    EulerDomainConstFieldView input{ 32, 5, values };

    EXPECT_NO_THROW( ValidateEulerDomainProblem( problem ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, input ) );
}

TEST( EulerDomainContract, RejectsInternalFieldShapeMismatch )
{
    const EulerDomainProblem problem = Problem();
    Real values[ 3 * 32 ] = {};
    EulerDomainConstFieldView wrongCells{ 31, 3, values };
    EulerDomainConstFieldView wrongEquations{ 32, 5, values };
    EulerDomainConstFieldView nullValues{ 32, 3, nullptr };

    EXPECT_THROW(
        ValidateEulerDomainField( problem, wrongCells ), std::invalid_argument );
    EXPECT_THROW(
        ValidateEulerDomainField( problem, wrongEquations ), std::invalid_argument );
    EXPECT_THROW(
        ValidateEulerDomainField( problem, nullValues ), std::invalid_argument );
}

TEST( EulerDomainContract, StateKeySeparatesExecutionOwnership )
{
    const EulerDomainStateKey cpu{ 0, 4, 0, AccelBackendKind::CPU, -1 };
    const EulerDomainStateKey same{ 0, 4, 0, AccelBackendKind::CPU, -1 };
    const EulerDomainStateKey otherBackend{ 0, 4, 0, AccelBackendKind::HIP, 0 };
    const EulerDomainStateKey otherDevice{ 0, 4, 0, AccelBackendKind::HIP, 1 };
    const EulerDomainStateKey otherZone{ 0, 5, 0, AccelBackendKind::CPU, -1 };

    EXPECT_TRUE( cpu == same );
    EXPECT_FALSE( cpu == otherBackend );
    EXPECT_FALSE( otherBackend == otherDevice );
    EXPECT_FALSE( cpu == otherZone );
}

TEST( AccelViewsContract, DescribesLayoutsAreaPolicyAndCapabilities )
{
    const SolverDomainCapabilities threeEq{
        3, 2, 1, true, true, true, true };
    const SolverDomainCapabilities fiveEq{
        5, 4, 2, true, true, true, true };
    const SolverDomainCapabilities fourEq{
        4, 0, 0, true, true, false, false };

    EXPECT_TRUE( threeEq.SupportsEulerState() );
    EXPECT_TRUE( fiveEq.SupportsEulerState() );
    EXPECT_FALSE( fourEq.SupportsEulerState() );

    Real values[ 3 * 4 ] = {};
    const SolverConstFieldView field{
        4, 3, values, FieldLayout::EquationMajor,
        FieldRepresentation::Primitive };
    EXPECT_NO_THROW( ValidateSolverFieldView( field ) );

    FaceGeometryView geometry;
    geometry.nFaces = 4;
    geometry.faceArea = values;
    geometry.areaPolicy = FaceAreaPolicy::BackendMultiplies;
    EXPECT_NO_THROW( ValidateFaceGeometryView( geometry ) );
}

TEST( AccelViewsContract, ValidatesCellGradientLayoutAndGhostExtent )
{
    Real q[ 4 ] = {};
    Real faceValues[ 3 ] = {};
    Real cellValues[ 4 ] = {};
    Real volume[ 2 ] = { 1.0, 1.0 };
    int left[ 3 ] = { 0, 1, 0 };
    int right[ 3 ] = { 2, 3, 1 };
    CellGradientView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nFaces = 3;
    view.nBoundaryFaces = 2;
    view.nEquations = 3;
    view.xFace = faceValues;
    view.yFace = faceValues;
    view.zFace = faceValues;
    view.xNormal = faceValues;
    view.yNormal = faceValues;
    view.zNormal = faceValues;
    view.faceArea = faceValues;
    view.xCell = cellValues;
    view.yCell = cellValues;
    view.zCell = cellValues;
    view.cellVolume = volume;
    view.leftCell = left;
    view.rightCell = right;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = q;
        view.dqdx[ equation ] = q;
        view.dqdy[ equation ] = q;
        view.dqdz[ equation ] = q;
    }
    EXPECT_NO_THROW( ValidateCellGradientView( view ) );
    view.nGhostCells = 1;
    EXPECT_THROW( ValidateCellGradientView( view ), std::invalid_argument );
}

TEST( AccelViewsContract, ValidatesDeviceReconstructionContract )
{
    Real cellValues[ 4 ] = {};
    Real faceValues[ 3 ] = {};
    int left[ 3 ] = { 0, 1, 0 };
    int right[ 3 ] = { 2, 3, 1 };
    unsigned char boundaryMask[ 3 ] = { 1, 1, 0 };
    ReconstructionBoundaryOperation boundaryOperation[ 3 ] = {
        ReconstructionBoundaryOperation::Preserve,
        ReconstructionBoundaryOperation::SolidOverride,
        ReconstructionBoundaryOperation::Preserve };
    BoundaryGhostUpdatePolicy boundaryGhostPolicy[ 3 ] = {
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::LocalBoundaryModel,
        BoundaryGhostUpdatePolicy::ExternalHalo };
    int owner = 0;

    CellFaceReconstructionView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nFaces = 3;
    view.nBoundaryFaces = 2;
    view.nEquations = 5;
    view.xFace = faceValues;
    view.yFace = faceValues;
    view.zFace = faceValues;
    view.xCell = cellValues;
    view.yCell = cellValues;
    view.zCell = cellValues;
    view.leftCell = left;
    view.rightCell = right;
    view.boundaryMask = boundaryMask;
    view.boundaryOperation = boundaryOperation;
    view.boundaryGhostPolicy = boundaryGhostPolicy;
    view.memorySpace = MemorySpace::Device;
    view.limiterMode = ReconstructionLimiterMode::Disabled;
    view.hasSolidBoundary = true;
    view.ownerToken = &owner;
    view.topologyGeneration = 2;
    view.fieldGeneration = 3;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = cellValues;
        view.dqdx[ equation ] = cellValues;
        view.dqdy[ equation ] = cellValues;
        view.dqdz[ equation ] = cellValues;
        view.bcQ[ equation ] = faceValues;
        view.qLeft[ equation ] = faceValues;
        view.qRight[ equation ] = faceValues;
    }

    EXPECT_NO_THROW( ValidateCellFaceReconstructionView( view ) );
    view.bcQ[ 4 ] = nullptr;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
}

TEST( AccelViewsContract, RejectsStaleOrIncompleteReconstructionContract )
{
    Real cellValues[ 4 ] = {};
    Real faceValues[ 3 ] = {};
    int leftCells[ 3 ] = { 0, 1, 0 };
    int rightCells[ 3 ] = { 2, 3, 1 };
    ReconstructionBoundaryOperation boundaryOperation[ 3 ] = {};
    BoundaryGhostUpdatePolicy boundaryGhostPolicy[ 3 ] = {
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::ExternalHalo };
    int owner = 0;
    CellFaceReconstructionView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nFaces = 3;
    view.nBoundaryFaces = 2;
    view.nEquations = 3;
    view.xFace = faceValues;
    view.yFace = faceValues;
    view.zFace = faceValues;
    view.xCell = cellValues;
    view.yCell = cellValues;
    view.zCell = cellValues;
    view.leftCell = leftCells;
    view.rightCell = rightCells;
    view.boundaryOperation = boundaryOperation;
    view.boundaryGhostPolicy = boundaryGhostPolicy;
    view.ownerToken = &owner;
    view.topologyGeneration = 1;
    view.fieldGeneration = 1;
    view.physicality = ReconstructionPhysicalityPolicy::None;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = cellValues;
        view.dqdx[ equation ] = cellValues;
        view.dqdy[ equation ] = cellValues;
        view.dqdz[ equation ] = cellValues;
        view.qLeft[ equation ] = faceValues;
        view.qRight[ equation ] = faceValues;
    }

    EXPECT_NO_THROW( ValidateCellFaceReconstructionView( view ) );
    view.topologyGeneration = 0;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
    view.topologyGeneration = 1;
    view.limiterMode = ReconstructionLimiterMode::Cell;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
}

TEST( AccelViewsContract, ValidatesReconstructionBoundaryOwnership )
{
    Real cellValues[ 4 ] = {};
    Real faceValues[ 3 ] = {};
    int left[ 3 ] = { 0, 1, 0 };
    int right[ 3 ] = { 2, 1, 3 };
    unsigned char boundaryMask[ 3 ] = { 1, 0, 1 };
    ReconstructionBoundaryOperation operations[ 3 ] = {};
    BoundaryGhostUpdatePolicy boundaryGhostPolicy[ 3 ] = {
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::ExternalHalo };
    int owner = 0;

    CellFaceReconstructionView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nFaces = 3;
    view.nBoundaryFaces = 2;
    view.nEquations = 3;
    view.xFace = faceValues;
    view.yFace = faceValues;
    view.zFace = faceValues;
    view.xCell = cellValues;
    view.yCell = cellValues;
    view.zCell = cellValues;
    view.leftCell = left;
    view.rightCell = right;
    view.boundaryMask = boundaryMask;
    view.boundaryOperation = operations;
    view.boundaryGhostPolicy = boundaryGhostPolicy;
    view.ownerToken = &owner;
    view.topologyGeneration = 1;
    view.fieldGeneration = 1;
    view.physicality = ReconstructionPhysicalityPolicy::None;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = cellValues;
        view.dqdx[ equation ] = cellValues;
        view.dqdy[ equation ] = cellValues;
        view.dqdz[ equation ] = cellValues;
        view.qLeft[ equation ] = faceValues;
        view.qRight[ equation ] = faceValues;
    }

    EXPECT_NO_THROW( ValidateCellFaceReconstructionView( view ) );

    right[ 0 ] = 1;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
    right[ 0 ] = 2;

    boundaryMask[ 1 ] = 1;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
    boundaryMask[ 1 ] = 0;

    left[ 1 ] = 2;
    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );
}

TEST( AccelViewsContract, RejectsAmbiguousGhostOwnership )
{
    Real cellValues[ 4 ] = {};
    Real faceValues[ 2 ] = {};
    int left[ 2 ] = { 0, 0 };
    int right[ 2 ] = { 1, 0 };
    ReconstructionBoundaryOperation operations[ 2 ] = {
        ReconstructionBoundaryOperation::Average,
        ReconstructionBoundaryOperation::Preserve };
    BoundaryGhostUpdatePolicy policies[ 2 ] = {
        BoundaryGhostUpdatePolicy::ExternalHalo,
        BoundaryGhostUpdatePolicy::ExternalHalo };
    int owner = 0;

    CellFaceReconstructionView view;
    view.nCells = 1;
    view.nGhostCells = 1;
    view.nFaces = 2;
    view.nBoundaryFaces = 1;
    view.nEquations = 3;
    view.xFace = faceValues;
    view.yFace = faceValues;
    view.zFace = faceValues;
    view.xCell = cellValues;
    view.yCell = cellValues;
    view.zCell = cellValues;
    view.leftCell = left;
    view.rightCell = right;
    view.boundaryOperation = operations;
    view.boundaryGhostPolicy = policies;
    view.ownerToken = &owner;
    view.topologyGeneration = 1;
    view.fieldGeneration = 1;
    view.physicality = ReconstructionPhysicalityPolicy::None;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = cellValues;
        view.dqdx[ equation ] = cellValues;
        view.dqdy[ equation ] = cellValues;
        view.dqdz[ equation ] = cellValues;
        view.qLeft[ equation ] = faceValues;
        view.qRight[ equation ] = faceValues;
    }

    EXPECT_THROW(
        ValidateCellFaceReconstructionView( view ), std::invalid_argument );

    policies[ 0 ] = BoundaryGhostUpdatePolicy::LocalBoundaryModel;
    operations[ 0 ] = ReconstructionBoundaryOperation::Average;
    EXPECT_NO_THROW( ValidateCellFaceReconstructionView( view ) );

    policies[ 0 ] = BoundaryGhostUpdatePolicy::ExternalHalo;
    operations[ 0 ] = ReconstructionBoundaryOperation::Preserve;
    EXPECT_NO_THROW( ValidateCellFaceReconstructionView( view ) );
}

TEST( AccelViewsContract, RejectsUnsupportedStateAndGeometryContracts )
{
    Real values[ 4 ] = {};
    int cells[ 4 ] = {};
    FaceStateView state;
    state.nFaces = 4;
    state.nEquations = 3;
    state.qLeft = values;
    state.qRight = values;
    state.layout = FieldLayout::EntityMajor;
    EXPECT_THROW( ValidateFaceStateView( state ), std::invalid_argument );

    FaceConnectivityView connectivity;
    connectivity.nFaces = 4;
    connectivity.nBoundaryFaces = 5;
    connectivity.leftCell = nullptr;
    connectivity.rightCell = cells;
    EXPECT_THROW( ValidateFaceConnectivityView( connectivity ),
        std::invalid_argument );

    FaceGeometryView geometry;
    geometry.nFaces = 4;
    EXPECT_THROW( ValidateFaceGeometryView( geometry ), std::invalid_argument );
}

TEST( AccelViewsContract, ValidatesCellStateUpdateContract )
{
    Real timeStep[ 2 ] = { 0.1, 0.2 };
    Real volume[ 2 ] = { 1.0, 1.0 };
    Real primitive[ 5 ][ 2 ] = {};
    int owner = 0;

    CellStateUpdateView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nEquations = 5;
    view.timeStep = timeStep;
    view.cellVolume = volume;
    view.gamma = 1.4;
    view.rkCoefficient = 0.5;
    view.cacheKey = &owner;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
        view.primitive[ equation ] = primitive[ equation ];

    EXPECT_NO_THROW( ValidateCellStateUpdateView( view ) );

    view.nEquations = 3;
    EXPECT_THROW(
        ValidateCellStateUpdateView( view ), std::invalid_argument );
    view.nEquations = 5;
    view.rkCoefficient = 0.0;
    EXPECT_THROW(
        ValidateCellStateUpdateView( view ), std::invalid_argument );
    view.rkCoefficient = 0.5;
    view.primitive[ 4 ] = nullptr;
    EXPECT_THROW(
        ValidateCellStateUpdateView( view ), std::invalid_argument );
}

TEST( AccelViewsContract, ValidatesCellTimeStepContract )
{
    Real primitive[ 5 ][ 4 ] = {};
    Real normal[ 3 ] = {};
    Real meshVelocity[ 3 ] = {};
    Real area[ 3 ] = { 1.0, 1.0, 1.0 };
    Real volume[ 2 ] = { 1.0, 1.0 };
    Real spectrum[ 2 ] = {};
    Real timeStep[ 2 ] = {};
    int leftCell[ 3 ] = { 0, 1, 0 };
    int rightCell[ 3 ] = { 2, 3, 1 };
    int owner = 0;

    CellTimeStepView view;
    view.nCells = 2;
    view.nGhostCells = 2;
    view.nFaces = 3;
    view.nBoundaryFaces = 2;
    view.nEquations = 5;
    view.xNormal = normal;
    view.yNormal = normal;
    view.zNormal = normal;
    view.meshVelocityNormal = meshVelocity;
    view.faceArea = area;
    view.cellVolume = volume;
    view.leftCell = leftCell;
    view.rightCell = rightCell;
    view.inviscidSpectrum = spectrum;
    view.timeStep = timeStep;
    view.gamma = 1.4;
    view.cfl = 0.01;
    view.cacheKey = &owner;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
        view.q[ equation ] = primitive[ equation ];

    EXPECT_NO_THROW( ValidateCellTimeStepView( view ) );

    view.cfl = 0.0;
    EXPECT_THROW(
        ValidateCellTimeStepView( view ), std::invalid_argument );
    view.cfl = 0.01;
    view.nGhostCells = 1;
    EXPECT_THROW(
        ValidateCellTimeStepView( view ), std::invalid_argument );
    view.nGhostCells = 2;
    view.q[ 4 ] = nullptr;
    EXPECT_THROW(
        ValidateCellTimeStepView( view ), std::invalid_argument );
}

} // namespace
