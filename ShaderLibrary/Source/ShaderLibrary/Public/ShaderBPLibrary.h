// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "RHICommandList.h"
#include "ShaderParameterStruct.h"
#include "Engine/TextureRenderTarget2D.h"
#include "RenderGraphUtils.h"
#include "ShaderBPLibrary.generated.h"

// UE4 compute shader tutorial
// https://github.com/Temaran/UE4ShaderPluginDemo/blob/4.25/Plugins/TemaranShaderTutorial/Source/ShaderDeclarationDemo/Private/ComputeShaderExample.cpp


class FPostProcessCS : public FGlobalShader
{
	DECLARE_GLOBAL_SHADER(FPostProcessCS);
	SHADER_USE_PARAMETER_STRUCT(FPostProcessCS, FGlobalShader);

	BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
		// You can only use InputTexture.Load() in shader if you do not provide a sampler parameter.
		SHADER_PARAMETER_SRV(Texture2D, InputTexture)
		SHADER_PARAMETER_UAV(RWTexture2D<float4>, OutputTexture)
		// Metal doesn't support GetDimensions(), so we send in this data via our parameters.
		SHADER_PARAMETER(FVector2D, TextureSize) 
	END_SHADER_PARAMETER_STRUCT()

public:

	//FPostProcessCS() {}

	//FPostProcessCS(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
	//	:FGlobalShader(Initializer) {
	//}

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters) {
		return FGlobalShader::ShouldCompilePermutation(Parameters);
	}

	#define NUM_THREADS_PER_GROUP_DIMENSION 32

	static void ModifyCompilationEnvironment(const FGlobalShaderPermutationParameters& Parameters, FShaderCompilerEnvironment& OutEnvironment) {
		FGlobalShader::ModifyCompilationEnvironment(Parameters, OutEnvironment);
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_X"), NUM_THREADS_PER_GROUP_DIMENSION);
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_Y"), NUM_THREADS_PER_GROUP_DIMENSION);
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_Z"), 1);
	}

};

BEGIN_SHADER_PARAMETER_STRUCT(FParticleSystemParameters, )
	
	// Simulation boundary
	SHADER_PARAMETER(FVector, SimCenter)
	SHADER_PARAMETER(FVector, SimHalfExtent)
	SHADER_PARAMETER(FVector, SimExtent)
	SHADER_PARAMETER(FVector, SimLowerCorner)
	SHADER_PARAMETER(FVector, SimHigherCorner)

	// Particle spawning
	SHADER_PARAMETER(FVector, SpawnCenter)
	SHADER_PARAMETER(FVector, SpawnHalfExtent)
	SHADER_PARAMETER(FVector, SpawnExtent)
	SHADER_PARAMETER(FVector, SpawnLowerCorner)
	SHADER_PARAMETER(FVector, SpawnHigherCorner)
	SHADER_PARAMETER(float, SpawnSpeed)
	SHADER_PARAMETER(int, ParticleNum)
	SHADER_PARAMETER(int, ActiveNum)
	SHADER_PARAMETER(int, FluidParticleNum)
	SHADER_PARAMETER(int, SolidParticleNum)

	// Physical quantities
	SHADER_PARAMETER(FVector, g)
	SHADER_PARAMETER(float, dt)
	SHADER_PARAMETER(float, rho)

	// Debugging
	SHADER_PARAMETER(float, SpeedScale)
	SHADER_PARAMETER(float, PressureScale)
	SHADER_PARAMETER(float, ViscosityScale)
	SHADER_PARAMETER(float, Gamma)
	SHADER_PARAMETER(float, Damping)
	SHADER_PARAMETER(int, ClampDensity)
	SHADER_PARAMETER(int, ClampPressure)
	SHADER_PARAMETER(int, ClampVelocity)
	SHADER_PARAMETER(int, ClampAcceleration)
	SHADER_PARAMETER(int, Preset)
	SHADER_PARAMETER(FVector, InitialVelocity)
	SHADER_PARAMETER(int, Visualization)
	SHADER_PARAMETER(float, Gamma1)
	SHADER_PARAMETER(float, Gamma2)
	SHADER_PARAMETER(int, UseNaiveBoundaryHandling)

	// SPH
	SHADER_PARAMETER(float, m)
	SHADER_PARAMETER(float, h)
	SHADER_PARAMETER(float, ParticleSpacing)
	SHADER_PARAMETER(int, ParticleNumPerLine)
	SHADER_PARAMETER(float, YoungsModulus)
	SHADER_PARAMETER(float, PoissonRatio)
	SHADER_PARAMETER(float, mu)
	SHADER_PARAMETER(float, lambda)
	SHADER_PARAMETER(float, SolidParticleMass)
	SHADER_PARAMETER(float, SolidParticleVolume)
	SHADER_PARAMETER(float, SolidParticleDensity)
	SHADER_PARAMETER(float, alpha)

	// Unit test
	SHADER_PARAMETER(int, UnitTestPreset)
	SHADER_PARAMETER(int, UnitTestNum)

	// Buffers
	// We use <float4> instead of <float3> since there is no DXGI_FORMAT_R32G32B32_FLOAT format in UE4
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4>, Position)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4>, Velocity)
	SHADER_PARAMETER_UAV(RWBuffer<float>, Mass)
	SHADER_PARAMETER_UAV(RWBuffer<float>, Volume)
	SHADER_PARAMETER_UAV(RWBuffer<float>, Density)
	SHADER_PARAMETER_UAV(RWBuffer<float>, Pressure)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4>, Acceleration)
	SHADER_PARAMETER_SRV(StructuredBuffer<float4>, RestPosition)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, J)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, L)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, RL)

	// Unit test
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, _A)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, _U)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, _Sigma)
	SHADER_PARAMETER_UAV(RWStructuredBuffer<float4x4>, _V)

END_SHADER_PARAMETER_STRUCT()


#define NUM_THREADS_PER_GROUP_DIMENSION 32

#define DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(ClassName) \
class ClassName : public FGlobalShader\
{\
	DECLARE_GLOBAL_SHADER(ClassName);\
	SHADER_USE_PARAMETER_STRUCT(ClassName, FGlobalShader);\
\
public:\
\
	using FParameters = FParticleSystemParameters;\
\
	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters) {\
		return FGlobalShader::ShouldCompilePermutation(Parameters);\
	}\
\
\
	static void ModifyCompilationEnvironment(const FGlobalShaderPermutationParameters& Parameters, FShaderCompilerEnvironment& OutEnvironment) {\
		FGlobalShader::ModifyCompilationEnvironment(Parameters, OutEnvironment);\
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_X"), NUM_THREADS_PER_GROUP_DIMENSION);\
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_Y"), 1);\
		OutEnvironment.SetDefine(TEXT("THREADGROUPSIZE_Z"), 1);\
	}\
\
};
// Can not put '\' at the last line of a macro

#define DECLARE_IMPLEMENT_DISPATCH_FUNCTION(ClassNameWithoutPrefix, ThreadNum)\
static void Dispatch##ClassNameWithoutPrefix##_RenderThread(FRHICommandListImmediate& RHICmdList, F##ClassNameWithoutPrefix::FParameters ShaderParameters) {\
	TShaderMapRef<F##ClassNameWithoutPrefix> ComputeShader(GetGlobalShaderMap(GMaxRHIFeatureLevel));\
	FIntVector GroupCount = FIntVector(FMath::DivideAndRoundUp(ThreadNum, NUM_THREADS_PER_GROUP_DIMENSION), 1, 1);\
	FComputeShaderUtils::Dispatch(RHICmdList, ComputeShader, ShaderParameters, GroupCount);\
};



DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FParticleSystemInitializationCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FParticleSystemSimulationCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHDensityCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHPressureForceCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHViscosityForceCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHAdvanceCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHSolidVolumeInitializationCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHSolidCorrectMatrixInitializationCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHSolidDeformationGradientCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FSPHSolidForceCS);
DECLARE_COMPUTE_SHADER_CLASS_USE_PARTICLE_SYSTEM_PARAMETERS(FUnitTestCS);


UCLASS(Blueprintable)
class SHADERLIBRARY_API UShaderBPLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()

public:
	//UFUNCTION() // FRHICommandListImmediate parameter can not be a UFUNCTION
	// The parameters of this function should not be FPostProcessCS::FParameters type.
	// You should use F-Class as parameters becuase that makes it convenient to manipulate data in RenderThread. 
	// Otherwise you need a more wrapper to manipulate render resources in RenderThread (Convert render resources to FPostProcessCS::FParameters).
	// (You can not manipulate render resources in GameThread, e.g. You can not use UTextureRenderTarget2D::GetRenderTargetResource()->GetRenderTargetUAV() in GameThread.)
	static void DispatchPostProcessCS_RenderThread(FRHICommandListImmediate& RHICmdList, UTexture* InputTexture, UTextureRenderTarget2D* OutputTexture);
	
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "DispatchPostProcessComputeShader"))
	static void DispatchPostProcessCS_GameThread(UTexture* InputTexture, UTextureRenderTarget2D* OutputTexture);

	// TODO: Chang RWBuffer<float4> to RWBuffer<float3>
	//static void DispatchParticleSystemInitializationCS_RenderThread(FRHICommandListImmediate& RHICmdList, int32 ParticleNum, FRWBuffer PositionBuffer);
	

	// No GameThread function version because I do not know how to pass FRWBuffer from GameThread to RenderThread
	//UFUNCTION(BlueprintCallable, meta = (DisplayName = "DispatchParticleSystemInitializationComputeShader"))
	//static void DispatchParticleSystemInitializationCS_GameThread();
	

	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(ParticleSystemInitializationCS, ShaderParameters.ParticleNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(ParticleSystemSimulationCS, ShaderParameters.ActiveNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHDensityCS, ShaderParameters.ActiveNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHPressureForceCS, ShaderParameters.ActiveNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHViscosityForceCS, ShaderParameters.ActiveNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHAdvanceCS, ShaderParameters.ActiveNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHSolidVolumeInitializationCS, ShaderParameters.SolidParticleNum - ShaderParameters.FluidParticleNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHSolidCorrectMatrixInitializationCS, ShaderParameters.SolidParticleNum - ShaderParameters.FluidParticleNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHSolidDeformationGradientCS, ShaderParameters.SolidParticleNum - ShaderParameters.FluidParticleNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(SPHSolidForceCS, ShaderParameters.SolidParticleNum - ShaderParameters.FluidParticleNum)
	DECLARE_IMPLEMENT_DISPATCH_FUNCTION(UnitTestCS, ShaderParameters.UnitTestNum)
	
};
