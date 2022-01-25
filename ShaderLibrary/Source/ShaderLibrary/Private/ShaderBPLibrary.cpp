// Fill out your copyright notice in the Description page of Project Settings.


#include "ShaderBPLibrary.h"

//							ClassName					ShaderPath				MainFunctionName	ShaderType
IMPLEMENT_GLOBAL_SHADER(FPostProcessCS,		"/ShaderLibrary/PostProcess.usf",		"MainCS",		SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FParticleSystemInitializationCS, "/ShaderLibrary/ParticleSystem.usf", "ParticleSystemIntializationCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FParticleSystemSimulationCS, "/ShaderLibrary/ParticleSystem.usf", "ParticleSystemSimulationCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHDensityCS, "/ShaderLibrary/ParticleSystem.usf", "SPHDensityCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHPressureForceCS, "/ShaderLibrary/ParticleSystem.usf", "SPHPressureForceCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHViscosityForceCS, "/ShaderLibrary/ParticleSystem.usf", "SPHViscosityForceCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHAdvanceCS, "/ShaderLibrary/ParticleSystem.usf", "SPHAdvanceCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHSolidVolumeInitializationCS, "/ShaderLibrary/ParticleSystem.usf", "SPHSolidVolumeInitializationCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHSolidCorrectMatrixInitializationCS, "/ShaderLibrary/ParticleSystem.usf", "SPHSolidCorrectMatrixInitializationCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHSolidDeformationGradientCS, "/ShaderLibrary/ParticleSystem.usf", "SPHSolidDeformationGradientCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FSPHSolidForceCS, "/ShaderLibrary/ParticleSystem.usf", "SPHSolidForceCS", SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FUnitTestCS, "/ShaderLibrary/ParticleSystem.usf", "UnitTestCS", SF_Compute);

void UShaderBPLibrary::DispatchPostProcessCS_RenderThread(FRHICommandListImmediate& RHICmdList, UTexture* InputTexture, UTextureRenderTarget2D* OutputTexture)
{
	check(IsInRenderingThread());

	// RDG can not use External UTextureRenderTarget2D, so, for simplicity, we do not use RDG.
	//FRDGBuilder GraphBuilder(RHICmdList);
	//GraphBuilder.RegisterExternalTexture(InputTexture);


	if (!OutputTexture->bCanCreateUAV) {
		// Set flag to permit creating UAVs for this resource, otherwise D3D can not create an UAV
		OutputTexture->bCanCreateUAV = true;
		// FTextureRenderTarget2DResource::InitDynamicRHI() will re-create resource using new bCanCreateUAV value.
		// New resource will be added to DeferredUpdateList, and will not be updated immediately.
		OutputTexture->GetRenderTargetResource()->UpdateRHI();
		// Update reource data immediately
		OutputTexture->UpdateResourceImmediate();
	}
	//FShaderResourceViewRHIRef InputTextureSRV = InputTexture->Resource->TextureRHI->GetNativeShaderResourceView();
	FRHITextureSRVCreateInfo SRVCreateInfo;
	FShaderResourceViewRHIRef InputTextureSRV = RHICreateShaderResourceView(InputTexture->Resource->TextureRHI, SRVCreateInfo);


	// Return nullptr
	//FUnorderedAccessViewRHIRef OutputTextureUAV = OutputTexture->GetRenderTargetResource()->GetRenderTargetUAV();

	// This Data is invalid (e.g. Texture Format is PF_UNKOWNN), and I do not know why
	//FUnorderedAccessViewRHIRef OutputTextureUAV = RHICreateUnorderedAccessView(OutputTexture->TextureReference.TextureReferenceRHI);
	
	// OK
	// There is no memory leak because of smart pointers.
	//FUnorderedAccessViewRHIRef OutputTextureUAV = RHICreateUnorderedAccessView(OutputTexture->GetRenderTargetResource()->GetRenderTargetTexture());
	
	// OK
	//FUnorderedAccessViewRHIRef OutputTextureUAV = RHICreateUnorderedAccessView(OutputTexture->GetRenderTargetResource()->TextureRHI);
	
	// OK
	FUnorderedAccessViewRHIRef OutputTextureUAV = RHICreateUnorderedAccessView(OutputTexture->Resource->TextureRHI);


	//// InputTexture is UTextureRenderTarget2D* type.
	//// We can find TexRef != (TexRef0 == TexRef1 == TexRef2 == TexRef3) by debugging their addresses ......
	//FTextureReferenceRHIRef TexRef = OutputTexture->TextureReference.TextureReferenceRHI.GetReference();
	//FTextureRHIRef TexRef0 = OutputTexture->TextureReference.TextureReferenceRHI.GetReference()->GetReferencedTexture();
	//FTextureRHIRef TexRef1 = OutputTexture->GetRenderTargetResource()->GetRenderTargetTexture();
	//FTextureRHIRef TexRef2 = OutputTexture->GetRenderTargetResource()->TextureRHI;
	//FTextureRHIRef TexRef3 = OutputTexture->Resource->TextureRHI;


	FPostProcessCS::FParameters ShaderParameters;
	ShaderParameters.InputTexture = InputTextureSRV;
	ShaderParameters.OutputTexture = OutputTextureUAV;
	ShaderParameters.TextureSize = FVector2D(InputTexture->GetSurfaceWidth(), OutputTexture->GetSurfaceHeight());


	TShaderMapRef<FPostProcessCS> ComputeShader(GetGlobalShaderMap(GMaxRHIFeatureLevel));

	FIntVector GroupCount = FIntVector(	FMath::DivideAndRoundUp(int(ShaderParameters.TextureSize.X), NUM_THREADS_PER_GROUP_DIMENSION), 
										FMath::DivideAndRoundUp(int(ShaderParameters.TextureSize.Y), NUM_THREADS_PER_GROUP_DIMENSION),
										1);

	FComputeShaderUtils::Dispatch(RHICmdList, ComputeShader, ShaderParameters, GroupCount);
}

void UShaderBPLibrary::DispatchPostProcessCS_GameThread(UTexture* InputTexture, UTextureRenderTarget2D* OutputTexture)
{
	ENQUEUE_RENDER_COMMAND(DispacthPostProcessCS)(
		[InputTexture, OutputTexture](FRHICommandListImmediate& RHICmdList)
		{
			DispatchPostProcessCS_RenderThread(RHICmdList, InputTexture, OutputTexture);
		});
}

//void UShaderBPLibrary::DispatchParticleSystemInitializationCS_RenderThread(FRHICommandListImmediate& RHICmdList, int32 ParticleNum, FRWBuffer PositionBuffer)
//{
//	FParticleSystemInitializationCS::FParameters ShaderParameters;
//	
//	// This format is PF_R32_FLOAT, this is not cosistent with the type (float4) in shader
//	//ShaderParameters.Position = PositionBuffer.UAV;
//
//	// This assignment will cause a crash. I don't know why
//	// ShaderParameters.Position = RHICreateUnorderedAccessView(PositionBuffer.Buffer, PF_A32B32G32R32F);;
//	
//	// We should use the correct format which is consistent with the type in shader
//	// RWBuffer<float4>  <--->   PF_A32B32G32R32F
//	// RWBuffer<float3>  <--->   (no such format in UE4)
//	FUnorderedAccessViewRHIRef PositionBufferUAV = RHICreateUnorderedAccessView(PositionBuffer.Buffer, PF_A32B32G32R32F);
//	ShaderParameters.Position = PositionBufferUAV;
//	
//	ShaderParameters.ParticleNum = ParticleNum;
//
//	TShaderMapRef<FParticleSystemInitializationCS> ComputeShader(GetGlobalShaderMap(GMaxRHIFeatureLevel));
//	FIntVector GroupCount = FIntVector(FMath::DivideAndRoundUp(ParticleNum, NUM_THREADS_PER_GROUP_DIMENSION), 1, 1);
//	FComputeShaderUtils::Dispatch(RHICmdList, ComputeShader, ShaderParameters, GroupCount);
//}

//void UShaderBPLibrary::DispatchParticleSystemInitializationCS_GameThread()
//{
//}
//
