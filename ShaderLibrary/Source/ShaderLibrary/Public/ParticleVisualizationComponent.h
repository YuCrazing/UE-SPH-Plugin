// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Components/InstancedStaticMeshComponent.h"
#include "ShaderBPLibrary.h"
#include "ParticleVisualizationComponent.generated.h"

UENUM()
enum class ESimulaitonType : uint8 {
	ParticleSystem,
	SPH,
	MPM
};

USTRUCT(Blueprintable)
struct FParticleSystemParameters_CPU{

	GENERATED_BODY()

	UPROPERTY(VisibleAnywhere)
	FVector SimCenter;
	UPROPERTY(VisibleAnywhere)
	FVector SimHalfExtent;
	UPROPERTY(VisibleAnywhere)
	FVector SimExtent;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	FVector SimLowerCorner;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	FVector SimHigherCorner;

	UPROPERTY(VisibleAnywhere)
	FVector SpawnCenter;
	FVector SpawnHalfExtent;
	UPROPERTY(VisibleAnywhere)
	FVector SpawnExtent;
	FVector SpawnLowerCorner;
	FVector SpawnHigherCorner;
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float SpawnSpeed;
	
	UPROPERTY(VisibleAnywhere)
	int32 ParticleNum;
	//UPROPERTY(EditAnywhere, BlueprintReadWrite)
	int32 ActiveNum;
	UPROPERTY(VisibleAnywhere)
	int32 FluidParticleNum;
	UPROPERTY(VisibleAnywhere)
	int32 SolidParticleNum;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	FVector g;
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float dt;
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float rho;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float SpeedScale;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float PressureScale;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float ViscosityScale;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float Gamma;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float Damping;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	uint8 ClampDensity : 1;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	uint8 ClampPressure : 1;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	uint8 ClampVelocity : 1;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	uint8 ClampAcceleration : 1;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	int Preset;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	FVector InitialVelocity;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	int32 Visualization;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float Gamma1;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Debug")
	float Gamma2;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	uint8 UseNaiveBoundaryHandling : 1;

	// SPH
	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float m;
	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float h;
	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float Volume;
	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float ParticleSpacing;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	int32 ParticleNumPerLine;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float ParticleRadius;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float ParticleVolumeScale;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float SmoothRadiusScale;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float YoungsModulus;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float PoissonRatio;

	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float mu;

	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float lambda;

	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float SolidParticleMass;

	UPROPERTY(VisibleAnywhere, Category = "SPH")
	float SolidParticleVolume;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float SolidParticleDensity;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SPH")
	float alpha;

	// Unit Test
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Unit Test")
	int32 UnitTestPreset;


	FParticleSystemParameters_CPU() {

		SimLowerCorner = FVector(-1, -1, -1);
		SimHigherCorner = FVector(1, 1, 1);

		SpawnSpeed = 1;
		
		g = FVector(0, 0, -9.8f);
		dt = 1.0f / 60;
		rho = 1000.0f;

		SpeedScale = 1.0f;
		PressureScale = 0.1f;
		ViscosityScale = 1.0f;
		Gamma = 1.0f;
		Damping = 0.99f;
		Preset = 3;
		Gamma1 = 0.7f;

		ParticleVolumeScale = 0.8f;
		SmoothRadiusScale = 4.0f;

		ReCalculateParameters();

	}

	void ReCalculateParameters() {

		ParticleSpacing = 2 * ParticleRadius;
		h = SmoothRadiusScale * ParticleRadius;

		float D = ParticleSpacing;

		float V = FMath::Pow(D, 3) * ParticleVolumeScale;

		m = V * rho;

		//ParticleNum = ParticleNumPerLine * ParticleNumPerLine * ParticleNumPerLine;
		

		SimCenter = (SimLowerCorner + SimHigherCorner) / 2;

		SimHalfExtent = SimHigherCorner - SimCenter;

		SimExtent = SimHalfExtent * 2;
		//SimLowerCorner = SimCenter - SimHalfExtent;
		//SimHigherCorner = SimCenter + SimHalfExtent;

		//ParticleNumPerLine = (int)FMath::Pow(ParticleNum + 0.5, 1.0 / 3);
		//SpawnExtent = FVector(FMath::Pow(Volume, 1.0 / 3));
		SpawnExtent = FVector(ParticleSpacing * (ParticleNumPerLine - 1));

		Volume = SpawnExtent.X * SpawnExtent.Y * SpawnExtent.Z;

		//SpawnCenter = SimCenter;


		//check(ParticleNumPerLine > 1);
		//ParticleSpacing = SpawnExtent.X / (ParticleNumPerLine - 1);

		//SpawnExtent = SpawnHalfExtent * 2;
		SpawnHalfExtent = SpawnExtent / 2;
		SpawnLowerCorner = SpawnCenter - SpawnHalfExtent;
		SpawnHigherCorner = SpawnCenter + SpawnHalfExtent;

		// Solid
		SolidParticleVolume = FMath::Pow(D, 3) * 0.8;
		SolidParticleMass = SolidParticleVolume * SolidParticleDensity;

		mu = YoungsModulus / (2 * (1 + PoissonRatio));
		lambda = YoungsModulus * PoissonRatio / ((1 + PoissonRatio) * (1 - 2 * PoissonRatio));
	}
};

USTRUCT(Blueprintable)
struct FFluidBlock {

	GENERATED_BODY()

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
		FVector LowerCorner;
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
		FVector HigherCorner;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
		FVector InitialVelocity;

	UPROPERTY(VisibleAnywhere, BlueprintReadWrite)
		FIntVector ParticleNumAlongAxis;
	UPROPERTY(VisibleAnywhere, BlueprintReadWrite)
		int ParticleNum;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
		uint8 Disable : 1;

	void ReCalculateParameters(const FParticleSystemParameters_CPU& Parameters) {
		ParticleNumAlongAxis = FIntVector((HigherCorner - LowerCorner) / Parameters.ParticleSpacing) + FIntVector(1);
		ParticleNum = ParticleNumAlongAxis.X * ParticleNumAlongAxis.Y * ParticleNumAlongAxis.Z;
	}
};


#define DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(ClassNameWithoutPrefix)\
void Dispatch##ClassNameWithoutPrefix##_GameThread()\
{\
	ENQUEUE_RENDER_COMMAND(ClassNameWithoutPrefix)(\
		[this](FRHICommandListImmediate& RHICmdList)\
	{\
		F##ClassNameWithoutPrefix::FParameters ShaderParameters; \
		CopyParticleSystemParameters(ShaderParameters); \
		UShaderBPLibrary::Dispatch##ClassNameWithoutPrefix##_RenderThread(RHICmdList, ShaderParameters); \
	}); \
};

UCLASS(Blueprintable, meta = (BlueprintSpawnableComponent))
class SHADERLIBRARY_API UParticleVisualizationComponent : public UInstancedStaticMeshComponent
{
	GENERATED_BODY()

	// We use FRWBuffer instead of FRWBufferStructured because we need a FVertexBufferRHIRef to replace this ISM's original vertex buffer.
	// Is it suitable putting FRWBuffer in GameThread?
	FRWBuffer PositionBuffer;
	FUnorderedAccessViewRHIRef PositionBufferUAV;
	TArray<FVector4> PositionBuffer_CPU;
	
	// TODO: change to FRWBufferStructured
	FRWBuffer VelocityBuffer;
	FUnorderedAccessViewRHIRef VelocityBufferUAV;
	TArray<FVector4> VelocityBuffer_CPU;

	FRWBuffer MassBuffer;

	FRWBuffer VolumeBuffer;

	FRWBuffer DensityBuffer;

	FRWBuffer PressureBuffer;

	// TODO: change to FRWBufferStructured
	FRWBuffer AccelerationBuffer;
	FUnorderedAccessViewRHIRef AccelerationBufferUAV;

	// TODO: change to FRWBufferStructured
	FRWBuffer RestPositionBuffer;
	FShaderResourceViewRHIRef RestPositionBufferSRV;

	FRWBufferStructured JBuffer;

	FRWBufferStructured LBuffer;
	FRWBufferStructured RLBuffer;

	FRWBufferStructured ABuffer;
	TArray<FMatrix> ABuffer_CPU;
	FRWBufferStructured UBuffer;
	FRWBufferStructured SigmaBuffer;
	FRWBufferStructured VBuffer;


public:

	UPROPERTY(VisibleAnywhere, BlueprintReadWrite, meta = (DisplayName = "Particle System Parameters"))
	FParticleSystemParameters_CPU Parameters;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (DisplayName = "Preset Parameters"))
	TArray<FParticleSystemParameters_CPU> Preset_Parameters;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	int32 SelectedParameterIndex;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	uint8 bUseCustomPositionBuffer : 1;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	uint8 bStartSimulate : 1;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	ESimulaitonType SimulationType;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	TArray<FFluidBlock> FluidBlock;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	TArray<FFluidBlock> SolidBlock;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	TArray<FFluidBlock> BoundaryBlock;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	int32 SubStep;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	uint8 RunUnitTest : 1;


	//UPROPERTY(VisibleAnywhere, BlueprintReadWrite)
	uint8 IsParameterValid : 1;

	float ActiveNum_Float;

	UParticleVisualizationComponent();

	void ReCalculateParameters();

	void CopyParticleSystemParameters(FParticleSystemInitializationCS::FParameters& ShaderParameters);

	void InitializeBuffers();

	void ReplaceInstanceData();

	virtual FPrimitiveSceneProxy* CreateSceneProxy() override;

	virtual void BeginPlay() override;
	virtual void TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;
	virtual void EndPlay(const EEndPlayReason::Type EndPlayReason) override;

#if WITH_EDITOR
	virtual void PostEditChangeChainProperty(FPropertyChangedChainEvent& PropertyChangedEvent) override;
#endif

	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(ParticleSystemInitializationCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(ParticleSystemSimulationCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHDensityCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHPressureForceCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHViscosityForceCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHAdvanceCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHSolidVolumeInitializationCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHSolidCorrectMatrixInitializationCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHSolidDeformationGradientCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(SPHSolidForceCS);
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION_GAMETHREAD(UnitTestCS);
	
};
