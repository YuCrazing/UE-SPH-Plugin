// Fill out your copyright notice in the Description page of Project Settings.


#include "ParticleVisualizationComponent.h"

#include "Engine/InstancedStaticMesh.h"
#include "ShaderBPLibrary.h"
#include "SceneManagement.h"

#define FLOAT_TYPE_SIZE_IN_BYTES 4
#define FLOATS_PER_POSITION 4
#define FLOATS_PER_VELOCITY 4
#define FLOATS_PER_ACCELERATION 4
#define FLOATS_PER_MATRIX_4X4 16

UParticleVisualizationComponent::UParticleVisualizationComponent()
{
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;
}

void UParticleVisualizationComponent::ReCalculateParameters()
{
	Parameters.ReCalculateParameters();
	Parameters.ParticleNum = 0;

	for (int32 i = 0; i < FluidBlock.Num(); ++i) {
		if (FluidBlock[i].Disable) continue;
		FluidBlock[i].ReCalculateParameters(Parameters);
		Parameters.ParticleNum += FluidBlock[i].ParticleNum;
	}
	Parameters.FluidParticleNum = Parameters.ParticleNum;

	for (int32 i = 0; i < SolidBlock.Num(); ++i) {
		if (SolidBlock[i].Disable) continue;
		SolidBlock[i].ReCalculateParameters(Parameters);
		Parameters.ParticleNum += SolidBlock[i].ParticleNum;
	}
	Parameters.SolidParticleNum = Parameters.ParticleNum;
	
	for (int32 i = 0; i < BoundaryBlock.Num(); ++i) {
		if (BoundaryBlock[i].Disable) continue;
		BoundaryBlock[i].ReCalculateParameters(Parameters);
		Parameters.ParticleNum += BoundaryBlock[i].ParticleNum;
	}

}

void UParticleVisualizationComponent::CopyParticleSystemParameters(FParticleSystemInitializationCS::FParameters& ShaderParameters)
{
	ShaderParameters.SimCenter = Parameters.SimCenter;
	ShaderParameters.SimHalfExtent = Parameters.SimHalfExtent;
	ShaderParameters.SimExtent = Parameters.SimExtent;
	ShaderParameters.SimLowerCorner = Parameters.SimLowerCorner;
	ShaderParameters.SimHigherCorner = Parameters.SimHigherCorner;

	ShaderParameters.SpawnCenter = Parameters.SpawnCenter;
	ShaderParameters.SpawnHalfExtent = Parameters.SpawnHalfExtent;
	ShaderParameters.SpawnExtent = Parameters.SpawnExtent;
	ShaderParameters.SpawnLowerCorner = Parameters.SpawnLowerCorner;
	ShaderParameters.SpawnHigherCorner = Parameters.SpawnHigherCorner;
	ShaderParameters.SpawnSpeed = Parameters.SpawnSpeed;
	ShaderParameters.ParticleNum = Parameters.ParticleNum;
	ShaderParameters.ActiveNum = Parameters.ActiveNum;
	ShaderParameters.FluidParticleNum = Parameters.FluidParticleNum;
	ShaderParameters.SolidParticleNum = Parameters.SolidParticleNum;

	ShaderParameters.g = Parameters.g;
	ShaderParameters.dt = Parameters.dt;
	ShaderParameters.rho = Parameters.rho;

	ShaderParameters.SpeedScale = Parameters.SpeedScale;
	ShaderParameters.PressureScale = Parameters.PressureScale;
	ShaderParameters.ViscosityScale = Parameters.ViscosityScale;
	ShaderParameters.Gamma = Parameters.Gamma;
	ShaderParameters.Damping = Parameters.Damping;
	ShaderParameters.ClampDensity = Parameters.ClampDensity;
	ShaderParameters.ClampPressure = Parameters.ClampPressure;
	ShaderParameters.ClampVelocity = Parameters.ClampVelocity;
	ShaderParameters.ClampAcceleration = Parameters.ClampAcceleration;
	ShaderParameters.Preset = Parameters.Preset;
	ShaderParameters.InitialVelocity = Parameters.InitialVelocity;
	ShaderParameters.Visualization = Parameters.Visualization;
	ShaderParameters.Gamma1 = Parameters.Gamma1;
	ShaderParameters.Gamma2 = Parameters.Gamma2;
	ShaderParameters.UseNaiveBoundaryHandling = Parameters.UseNaiveBoundaryHandling;

	ShaderParameters.m = Parameters.m;
	ShaderParameters.h = Parameters.h;
	ShaderParameters.ParticleSpacing = Parameters.ParticleSpacing;
	ShaderParameters.ParticleNumPerLine = Parameters.ParticleNumPerLine;
	ShaderParameters.YoungsModulus = Parameters.YoungsModulus;
	ShaderParameters.PoissonRatio = Parameters.PoissonRatio;
	ShaderParameters.mu = Parameters.mu;
	ShaderParameters.lambda = Parameters.lambda;
	ShaderParameters.SolidParticleMass = Parameters.SolidParticleMass;
	ShaderParameters.SolidParticleVolume = Parameters.SolidParticleVolume;
	ShaderParameters.SolidParticleDensity = Parameters.SolidParticleDensity;
	ShaderParameters.alpha = Parameters.alpha;

	// Unit test
	ShaderParameters.UnitTestPreset = Parameters.UnitTestPreset;
	ShaderParameters.UnitTestNum = ABuffer_CPU.Num();

	// Buffers

	// Using a temporary variable causes a crash because of null pointers
	//FUnorderedAccessViewRHIRef PositionBufferUAV = RHICreateUnorderedAccessView(PositionBuffer.Buffer, PF_A32B32G32R32F);

	ShaderParameters.Position = PositionBufferUAV;
	ShaderParameters.Velocity = VelocityBufferUAV;
	ShaderParameters.Mass = MassBuffer.UAV;
	ShaderParameters.Volume = VolumeBuffer.UAV;
	ShaderParameters.Density = DensityBuffer.UAV;
	ShaderParameters.Pressure = PressureBuffer.UAV;
	ShaderParameters.Acceleration = AccelerationBufferUAV;
	ShaderParameters.RestPosition = RestPositionBufferSRV;

	ShaderParameters.J = JBuffer.UAV;
	ShaderParameters.L = LBuffer.UAV;
	ShaderParameters.RL = RLBuffer.UAV;

	// Unit test
	ShaderParameters._A = ABuffer.UAV;
	ShaderParameters._U = UBuffer.UAV;
	ShaderParameters._Sigma = SigmaBuffer.UAV;
	ShaderParameters._V = VBuffer.UAV;

}

void UParticleVisualizationComponent::InitializeBuffers()
{
	ENQUEUE_RENDER_COMMAND(CustomPositionBuffer_Initialization)(
		[this](FRHICommandListImmediate& RHICmdList)
		{
			// We must make sure ParticleNum > 0 otherwise a crash will happen.
			PositionBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, FLOATS_PER_POSITION * Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("PositionBuffer"));
			PositionBufferUAV = RHICreateUnorderedAccessView(PositionBuffer.Buffer, PF_A32B32G32R32F);
		
			VelocityBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, FLOATS_PER_VELOCITY * Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("VelocityBuffer"));
			VelocityBufferUAV = RHICreateUnorderedAccessView(VelocityBuffer.Buffer, PF_A32B32G32R32F);

			MassBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("MassBuffer"));
			
			VolumeBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("VolumeBuffer"));
			
			DensityBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("DensityBuffer"));
			
			PressureBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("PressureBuffer"));
			
			AccelerationBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, FLOATS_PER_ACCELERATION * Parameters.ParticleNum, PF_R32_FLOAT, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("AccelerationBuffer"));
			AccelerationBufferUAV = RHICreateUnorderedAccessView(AccelerationBuffer.Buffer, PF_A32B32G32R32F);

			RestPositionBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES, FLOATS_PER_POSITION * Parameters.ParticleNum, PF_R32_FLOAT, BUF_ShaderResource, TEXT("RestPositionBuffer"));
			RestPositionBufferSRV = RHICreateShaderResourceView(RestPositionBuffer.Buffer, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_POSITION, PF_A32B32G32R32F);

			JBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4, Parameters.ParticleNum, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("JBuffer"));
			LBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4, Parameters.ParticleNum, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("LBuffer"));
			RLBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4, Parameters.ParticleNum, BUF_UnorderedAccess | BUF_ShaderResource, TEXT("RLBuffer"));


			ReplaceInstanceData();


			PositionBuffer_CPU.Empty();
			PositionBuffer_CPU.AddDefaulted(Parameters.ParticleNum);

			VelocityBuffer_CPU.Empty();
			VelocityBuffer_CPU.AddDefaulted(Parameters.ParticleNum);

			int32 IndexOffset = 0;
			for (int32 BlockIndex = 0; BlockIndex < FluidBlock.Num(); ++BlockIndex) {
				if (FluidBlock[BlockIndex].Disable) continue;
				FIntVector N = FluidBlock[BlockIndex].ParticleNumAlongAxis;
				for (int32 _i = 0; _i < FluidBlock[BlockIndex].ParticleNum; ++_i) {
					int32 i = IndexOffset + _i;
					PositionBuffer_CPU[i] = FVector4(FluidBlock[BlockIndex].LowerCorner + Parameters.ParticleSpacing * FVector(_i % N.X, _i / N.X % N.Y, _i / N.X / N.Y), 0);
					VelocityBuffer_CPU[i] = FVector4(FluidBlock[BlockIndex].InitialVelocity, 0);
				}
				IndexOffset += FluidBlock[BlockIndex].ParticleNum;
			}

			for (int32 BlockIndex = 0; BlockIndex < SolidBlock.Num(); ++BlockIndex) {
				if (SolidBlock[BlockIndex].Disable) continue;
				FIntVector N = SolidBlock[BlockIndex].ParticleNumAlongAxis;
				for (int32 _i = 0; _i < SolidBlock[BlockIndex].ParticleNum; ++_i) {
					int32 i = IndexOffset + _i;
					PositionBuffer_CPU[i] = FVector4(SolidBlock[BlockIndex].LowerCorner + Parameters.ParticleSpacing * FVector(_i % N.X, _i / N.X % N.Y, _i / N.X / N.Y), 0);
					VelocityBuffer_CPU[i] = FVector4(SolidBlock[BlockIndex].InitialVelocity, 0);
				}
				IndexOffset += SolidBlock[BlockIndex].ParticleNum;
			}

			for (int32 BlockIndex = 0; BlockIndex < BoundaryBlock.Num(); ++BlockIndex) {
				if (BoundaryBlock[BlockIndex].Disable) continue;
				FIntVector N = BoundaryBlock[BlockIndex].ParticleNumAlongAxis;
				for (int32 _i = 0; _i < BoundaryBlock[BlockIndex].ParticleNum; ++_i) {
					int32 i = IndexOffset + _i;
					PositionBuffer_CPU[i] = FVector4(BoundaryBlock[BlockIndex].LowerCorner + Parameters.ParticleSpacing * FVector(_i % N.X, _i / N.X % N.Y, _i / N.X / N.Y), 0);
					VelocityBuffer_CPU[i] = FVector4(BoundaryBlock[BlockIndex].InitialVelocity, 0);
				}
				IndexOffset += BoundaryBlock[BlockIndex].ParticleNum;
			}

			void* VertexBufferData = RHILockVertexBuffer(PositionBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_POSITION * Parameters.ParticleNum, RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, PositionBuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_POSITION * Parameters.ParticleNum);
			RHIUnlockVertexBuffer(PositionBuffer.Buffer);

			VertexBufferData = RHILockVertexBuffer(RestPositionBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_POSITION * Parameters.ParticleNum, RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, PositionBuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_POSITION * Parameters.ParticleNum);
			RHIUnlockVertexBuffer(RestPositionBuffer.Buffer);

			VertexBufferData = RHILockVertexBuffer(VelocityBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_VELOCITY * Parameters.ParticleNum, RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, VelocityBuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_VELOCITY * Parameters.ParticleNum);
			RHIUnlockVertexBuffer(VelocityBuffer.Buffer);


			// Unit test data


			ABuffer_CPU.Add(
				FMatrix(
					FPlane(0.0, 0.0, 0.0, 0.0),
					FPlane(0.0, 0.0, 0.0, 0.0),
					FPlane(0.0, 0.0, 0.0, 0.0),
					FPlane(0.0, 0.0, 0.0, 0.0)
				)
			);

			ABuffer_CPU.Add(
				FMatrix(
					FPlane(1.0, 0.0, 0.0, 0.0),
					FPlane(0.0, 1.0, 0.0, 0.0),
					FPlane(0.0, 0.0, 1.0, 0.0),
					FPlane(0.0, 0.0, 0.0, 0.0)
				)
			);
			ABuffer_CPU.Add(
				FMatrix(
					FPlane(8.0, -6.0, 2.0, 0.0),
					FPlane(-6.0, 7.0, -4.0, 0.0),
					FPlane(2.0, -4.0, 3.0, 0.0),
					FPlane(0.0, 0.0, 0.0, 0.0)
				)
			);
			/* answers:
			
				U:
				0.66666667	-0.66666667	0.33333333
				-0.66666667	-0.33333333	0.66666667
				0.33333333	0.66666667	0.66666667

				Sigma:
				15	0	0
				0	3	0
				0	0	0

				V:
				0.66666667	-0.66666667	0.33333333
				-0.66666667	-0.33333333	0.66666667
				0.33333333	0.66666667	0.66666667

				or

				U:
				-0.66666667	-0.66666667	0.33333333
				0.66666667	-0.33333333	0.66666667
				-0.33333333	0.66666667	0.66666667

				Sigma:
				15	0	0
				0	3	0
				0	0	0

				V:
				-0.66666667	-0.66666667	0.33333333
				0.66666667	-0.33333333	0.66666667
				-0.33333333	0.66666667	0.66666667
			*/



			ABuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES* FLOATS_PER_MATRIX_4X4, ABuffer_CPU.Num(), BUF_UnorderedAccess | BUF_ShaderResource, TEXT("ABuffer"));
			UBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES* FLOATS_PER_MATRIX_4X4, ABuffer_CPU.Num(), BUF_UnorderedAccess | BUF_ShaderResource, TEXT("UBuffer"));
			SigmaBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES* FLOATS_PER_MATRIX_4X4, ABuffer_CPU.Num(), BUF_UnorderedAccess | BUF_ShaderResource, TEXT("SigmaBuffer"));
			VBuffer.Initialize(FLOAT_TYPE_SIZE_IN_BYTES* FLOATS_PER_MATRIX_4X4, ABuffer_CPU.Num(), BUF_UnorderedAccess | BUF_ShaderResource, TEXT("VBuffer"));

			VertexBufferData = RHILockStructuredBuffer(ABuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num(), RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, ABuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num());
			RHIUnlockStructuredBuffer(ABuffer.Buffer);

			VertexBufferData = RHILockStructuredBuffer(UBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num(), RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, ABuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num());
			RHIUnlockStructuredBuffer(UBuffer.Buffer);

			VertexBufferData = RHILockStructuredBuffer(SigmaBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num(), RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, ABuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num());
			RHIUnlockStructuredBuffer(SigmaBuffer.Buffer);

			VertexBufferData = RHILockStructuredBuffer(VBuffer.Buffer, 0, FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num(), RLM_WriteOnly);
			FMemory::Memcpy(VertexBufferData, ABuffer_CPU.GetData(), FLOAT_TYPE_SIZE_IN_BYTES * FLOATS_PER_MATRIX_4X4 * ABuffer_CPU.Num());
			RHIUnlockStructuredBuffer(VBuffer.Buffer);


		});
}

void UParticleVisualizationComponent::ReplaceInstanceData()
{
	FVertexBufferRHIRef OldVertexBuffer = PerInstanceRenderData->InstanceBuffer.InstanceCustomDataBuffer.VertexBufferRHI;
	FVertexBufferRHIRef NewVertexBuffer = PositionBuffer.Buffer;

	if (OldVertexBuffer != NewVertexBuffer) {
		PerInstanceRenderData->InstanceBuffer.InstanceCustomDataBuffer.VertexBufferRHI = PositionBuffer.Buffer;
		PerInstanceRenderData->InstanceBuffer.InstanceCustomDataSRV = PositionBuffer.SRV;
		//InstanceUpdateCmdBuffer.Edit();
		MarkRenderStateDirty();
	}

}

FPrimitiveSceneProxy* UParticleVisualizationComponent::CreateSceneProxy()
{
	FPrimitiveSceneProxy* _SceneProxy = Super::CreateSceneProxy();
	return _SceneProxy;
}

void UParticleVisualizationComponent::BeginPlay()
{
	Super::BeginPlay();

	ActiveNum_Float = 0;

	ReCalculateParameters();


	IsParameterValid = (Parameters.ParticleNum > 0);
	
	if (!ensure(IsParameterValid)) {
		return;
	}


	// Set instance number
	ClearInstances();
	TArray<FTransform> ParticleTransforms;
	ParticleTransforms.AddDefaulted(Parameters.ParticleNum);
	AddInstances(ParticleTransforms, false);


	// We use <float4> instead of <float3> since there is no DXGI_FORMAT_R32G32B32_FLOAT format in UE4
	NumCustomDataFloats = FLOATS_PER_POSITION;
	// Clear out and reinit to 0
	PerInstanceSMCustomData.Empty(PerInstanceSMData.Num() * NumCustomDataFloats);
	PerInstanceSMCustomData.SetNumZeroed(PerInstanceSMData.Num() * NumCustomDataFloats);

	MarkRenderStateDirty();




	if (bUseCustomPositionBuffer) {
		InitializeBuffers();
		DispatchParticleSystemInitializationCS_GameThread();
		DispatchSPHSolidVolumeInitializationCS_GameThread();
		DispatchSPHSolidCorrectMatrixInitializationCS_GameThread();
	}
	else {
		for (int i = 0; i < Parameters.ParticleNum; ++i) {
			for (int j = 0; j < NumCustomDataFloats; ++j) {
				SetCustomDataValue(i, j, j, false);
			}
		}
		SetCustomDataValue(0, 0, 0, true);
	}



}

void UParticleVisualizationComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	if (!IsParameterValid) {
		return;
	}

	if (bUseCustomPositionBuffer) {

		if (bStartSimulate) {
			ActiveNum_Float += Parameters.SpawnSpeed;
			ActiveNum_Float = FMath::Min((int)ActiveNum_Float, (int)Parameters.ParticleNum);
			Parameters.ActiveNum = (int32)ActiveNum_Float;
			//Parameters.ActiveNum = -(int32)2834346832;
			Parameters.ActiveNum = FMath::Min(Parameters.ActiveNum, Parameters.ParticleNum);
			//TArray<FTransform> ParticleTransforms;
			//ParticleTransforms.AddDefaulted(Parameters.ParticleNum);
			//AddInstances(ParticleTransforms, false);

			for (int32 i = 0; i < SubStep; ++i) {
				switch (SimulationType)
				{
				case ESimulaitonType::ParticleSystem:
					DispatchParticleSystemSimulationCS_GameThread();
					break;
				case ESimulaitonType::SPH:
					DispatchSPHDensityCS_GameThread();
					DispatchSPHPressureForceCS_GameThread();
					DispatchSPHViscosityForceCS_GameThread();
					//DispatchSPHSolidVolumeInitializationCS_GameThread();
					//DispatchSPHSolidCorrectMatrixInitializationCS_GameThread();
					DispatchSPHSolidDeformationGradientCS_GameThread();
					DispatchSPHSolidForceCS_GameThread();
					DispatchSPHAdvanceCS_GameThread();
					break;
				case ESimulaitonType::MPM:
					break;
				default:
					break;
				}
			}
		}

		if (RunUnitTest) {
			DispatchUnitTestCS_GameThread();
		}

		ReplaceInstanceData();
	}
}

void UParticleVisualizationComponent::EndPlay(const EEndPlayReason::Type EndPlayReason)
{
	Super::EndPlay(EndPlayReason);

	ENQUEUE_RENDER_COMMAND(CustomVertexBuffer_Release)(
		[this](FRHICommandListImmediate& RHICmdList)
		{
			this->PositionBuffer.Release();
		});
}

#if WITH_EDITOR
void UParticleVisualizationComponent::PostEditChangeChainProperty(FPropertyChangedChainEvent& PropertyChangedEvent)
{
	if (SelectedParameterIndex >= 0 && SelectedParameterIndex < Preset_Parameters.Num()) {
		Parameters = Preset_Parameters[SelectedParameterIndex];
	}

	ReCalculateParameters();

	//Parameters.ReCalculateParameters();

	Super::PostEditChangeChainProperty(PropertyChangedEvent);
}
#endif
