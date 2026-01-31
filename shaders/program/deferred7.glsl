varying vec2 texcoord;

varying vec3 sunWorldDir, moonWorldDir, lightWorldDir;
varying vec3 sunViewDir, moonViewDir, lightViewDir;

varying float isNoon, isNight, sunRiseSet;
varying float isNoonS, isNightS, sunRiseSetS;

varying vec3 sunColor, skyColor;


#include "/lib/uniform.glsl"
#include "/lib/settings.glsl"
#include "/lib/common/utils.glsl"
#include "/lib/camera/colorToolkit.glsl"
#include "/lib/camera/filter.glsl"
#include "/lib/common/position.glsl"
#include "/lib/common/normal.glsl"
#include "/lib/common/noise.glsl"

#include "/lib/atmosphere/atmosphericScattering.glsl"

#ifdef FSH

const bool shadowtex0Mipmap = false;
const bool shadowtex1Mipmap = false;
const bool shadowcolor0Mipmap = false;
const bool shadowcolor1Mipmap = false;



#include "/lib/common/gbufferData.glsl"
#include "/lib/common/materialIdMapper.glsl"
#include "/lib/lighting/lightmap.glsl"
#include "/lib/lighting/shadowMapping.glsl"
#include "/lib/lighting/screenSpaceShadow.glsl"
#include "/lib/lighting/RSM.glsl"
#include "/lib/lighting/SSAO.glsl"
#include "/lib/lighting/localLightShadow.glsl"
#include "/lib/surface/PBR.glsl"

#include "/lib/atmosphere/fog.glsl"
#include "/lib/atmosphere/celestial.glsl"
// #include "/lib/atmosphere/volumetricClouds.glsl"

#ifdef VOLUMETRIC_CLOUDS
vec4 sampleCloudUpscaled(vec2 uv){
	vec4 c0 = texture(colortex3, uv);
	vec2 texel = invViewSize;
	vec4 c1 = texture(colortex3, uv + vec2(texel.x, 0.0));
	vec4 c2 = texture(colortex3, uv + vec2(0.0, texel.y));
	vec4 c3 = texture(colortex3, uv + texel);

	float edge = CLOUD_UPSAMPLE_EDGE;
	float w1 = 1.0 - saturate(abs(c1.a - c0.a) * edge);
	float w2 = 1.0 - saturate(abs(c2.a - c0.a) * edge);
	float w3 = 1.0 - saturate(abs(c3.a - c0.a) * edge);

	vec4 sum = c0 + c1 * w1 + c2 * w2 + c3 * w3;
	float wsum = 1.0 + w1 + w2 + w3;
	return sum / max(wsum, 1e-4);
}
#endif



void main() {
	vec4 CT2 = texelFetch(colortex2, ivec2(gl_FragCoord.xy), 0);
	vec4 color = texture(colortex0, texcoord);	// albedo
	vec3 texColor = color.rgb;
	vec3 albedo = pow(texColor, vec3(2.2));
	vec3 diffuse = albedo / PI;

	vec3 normalV = normalize(normalDecode(normalEnc));
	vec3 normalW = normalize(viewPosToWorldPos(vec4(normalV, 0.0)).xyz);

	vec3 L2 = BLACK;
	vec3 ao = vec3(1.0);
	vec3 tcltDbgRad = vec3(0.0);
	float tcltDbgConf = 0.0;
	vec3 tcltDbgColor = vec3(0.0);

	#if defined DISTANT_HORIZONS && !defined NETHER && !defined END
		bool isTerrain = skyB < 0.5;

		float depth1;
		vec4 viewPos1;
		if(dhTerrain > 0.5){ 
			float dhDepth = texture(dhDepthTex0, texcoord).r;
			viewPos1 = screenPosToViewPosDH(vec4(unTAAJitter(texcoord), dhDepth, 1.0));
			depth1 = viewPosToScreenPos(viewPos1).z;
		}else{
			depth1 = texture(depthtex1, texcoord).r;
			viewPos1 = screenPosToViewPos(vec4(unTAAJitter(texcoord), depth1, 1.0));	
		}
	#else 
		bool isTerrain = skyB < 0.5;

		float depth1 = texture(depthtex1, texcoord).r;
		vec4 viewPos1 = screenPosToViewPos(vec4(unTAAJitter(texcoord), depth1, 1.0));	
	#endif

	vec3 viewDir = normalize(viewPos1.xyz);
	vec4 worldPos1 = viewPosToWorldPos(viewPos1);
	vec3 worldDir = normalize(worldPos1.xyz);
	vec3 shadowPos = getShadowPos(worldPos1).xyz;
	float worldDis1 = length(worldPos1);

	vec4 viewPos1R = screenPosToViewPos(vec4(texcoord, depth1, 1.0));
	vec4 worldPos1R = viewPosToWorldPos(viewPos1R);
	vec2 prePos = getPrePos(worldPos1R).xy;
	vec2 velocity = texcoord - prePos;

	if(isTerrain){	
		vec2 lightmap = AdjustLightmap(mcLightmap);
		#ifdef HELD_BLOCK_DYNAMIC_LIGHT
			float heldBlockLight = max(heldBlockLightValue, heldBlockLightValue2) / 15.0;
			heldBlockLight *= pow(remapSaturate(worldDis1, 0.0, DYNAMIC_LIGHT_DISTANCE, 1.0, 0.0), ARTIFICIAL_LIGHT_FALLOFF);
			#ifdef HELD_BLOCK_NORMAL_AFFECT
				heldBlockLight *= saturate(dot(normalV, -normalize(vec3(viewPos1.xyz))));
			#endif
			lightmap.x = max(lightmap.x, heldBlockLight);
		#endif

		MaterialParams materialParams = MapMaterialParams(specularMap);
		#ifdef PBR_REFLECTIVITY
			mat2x3 PBR = CalculatePBR(viewDir, normalV, lightViewDir, albedo, materialParams);
			vec3 BRDF = PBR[0] + PBR[1];
			vec3 BRDF_D = BRDF_Diffuse(normalV, viewDir, albedo, materialParams);
		#else
			vec3 BRDF = albedo / PI;
			vec3 BRDF_D = BRDF;
		#endif



		float cos_theta_O = dot(normalW, lightWorldDir);
		float cos_theta = max(cos_theta_O, 0.0);

		// bzyzhang: 练习项目(十一)：次表面散射的近似实现
		// https://zhuanlan.zhihu.com/p/348106844
		float sssWrap = SSS_INTENSITY * materialParams.subsurfaceScattering;
		if(plants > 0.5) sssWrap = 20.0;
		cos_theta = saturate((cos_theta_O + sssWrap) / (1 + sssWrap));


		float noRSM = hand > 0.5 ? 1.0 : 0.0;
		float lightMask = pow(smoothstep(0.0, 0.1, lightmap.y), 0.45);
		float UoN = dot(normalW, upWorldDir);
		vec3 skyLight = lightmap.y * BRDF_D
					* mix(sunColor, skyColor, SUN_SKY_BLEND - 0.05 * noRSM * lightmap.y)
					* mix(1.0, UoN * 0.5 + 0.5, 0.75);
		

		vec4 gi = getGI(depth1, normalW);
		gi.a = 1.0 - gi.a;
		if(noRSM < 0.5) {
			L2 = sunColor * BRDF_D * gi.rgb;
			#ifdef AO_ENABLED
				#ifdef AO_MULTI_BOUNCE
					ao = AOMultiBounce(albedo, saturate(gi.a));
				#else 
					ao = vec3(saturate(gi.a));
				#endif
			#endif
		}



		float shadow = 1.0;
		vec3 colorShadow = vec3(0.0);
		if(!outScreen(shadowPos.xy) && cos_theta > 0.001){
			shadow = CT4R.x;
			colorShadow = getColorShadow(shadowPos, shadow);
		}
		float RTShadow = CT4R.y;
		shadow = min(shadow, RTShadow);

		shadow = saturate(shadow);
		vec3 visibility = vec3(shadow + colorShadow * 1.0);
		vec3 direct = sunColor * BRDF * visibility * cos_theta;



		vec3 artificialColor = artificial_color;
		#ifdef TCLT_COLORED_LIGHTING
			vec4 tclt = texture(colortex11, texcoord);
			vec3 tcltRad = tclt.rgb * TCLT_STRENGTH;
			float tcltConf = tclt.a;
			vec3 tcltLocalRad = vec3(0.0);
			float tcltLocalConf = 0.0;

			if (tcltConf < 0.25) {
				vec2 jitter = hash42(gl_FragCoord.xy).rg * 2.0 - 1.0;
				vec2 step = invViewSize * mix(8.0, 48.0, saturate(lightmap.x * 1.2));

				vec2 offsets[12];
				offsets[0]  = vec2( 0.75,  0.25);
				offsets[1]  = vec2(-0.25,  0.75);
				offsets[2]  = vec2( 0.25, -0.75);
				offsets[3]  = vec2(-0.75, -0.25);
				offsets[4]  = vec2( 1.25,  0.0);
				offsets[5]  = vec2(-1.25,  0.0);
				offsets[6]  = vec2( 0.0,  1.25);
				offsets[7]  = vec2( 0.0, -1.25);
				offsets[8]  = vec2( 1.75,  0.75);
				offsets[9]  = vec2(-1.75,  0.75);
				offsets[10] = vec2( 0.75, -1.75);
				offsets[11] = vec2(-0.75, -1.75);

				vec3 bestColor = vec3(0.0);
				float bestW = 0.0;
				float depthLin = linearizeDepth(depth1);
				const float blockIDRange = 0.3;

				for (int i = 0; i < 12; ++i) {
					vec2 uv = texcoord + (offsets[i] + jitter * 0.35) * step;
					if (outScreen(uv)) continue;

					float d = texture(depthtex1, uv).r;
					if (d >= 1.0) continue;
					float dLin = linearizeDepth(d);
					if (abs(dLin - depthLin) > TCLT_DEPTH_TH) continue;

					vec4 ct4 = texture(colortex4, uv);
					vec2 ct4g = unpack16To2x8(ct4.g);
					float bid = ct4g.x * ID_SCALE;
					float glow = checkInRange(bid, GLOWING_BLOCK, blockIDRange) + checkInRange(bid, NO_ANISO, blockIDRange);

					vec4 specS = unpack2x16To4x8(ct4.ba);
					float em = specS.a;
					if (em * 255.0 > 254.1) em = 0.0;

					vec3 al = pow(texture(colortex0, uv).rgb, vec3(2.2));
					vec2 lmS = texture(colortex5, uv).ba;
					float bl = lmS.x;
					float lum = getLuminance(al);
					float chroma = max3(al.r, al.g, al.b) - min3(al.r, al.g, al.b);
					float lumW = smoothstep(0.05, 0.6, lum);
					float chromaW = smoothstep(0.05, 0.35, chroma);

					float src = max(em * 2.0, glow * 1.5);
					float chromaSrc = bl * chromaW * 0.5;
					float weight = max(src, chromaSrc) * smoothstep(0.05, 0.8, bl);
					if (weight <= 0.01) continue;

					if (weight > bestW) {
						bestW = weight;
						bestColor = al;
					}
				}

				vec3 localRad = bestColor * bestW * TCLT_EMISSIVE_GAIN;
				tcltLocalRad = localRad;
				tcltLocalConf = saturate(bestW * 2.5);
				tcltRad = mix(localRad, tcltRad, tcltConf) * TCLT_STRENGTH;
				tcltConf = max(tcltConf, tcltLocalConf);
			}
			float tcltMax = max3(tcltRad.r, tcltRad.g, tcltRad.b);
			vec3 tcltColor = tcltMax > 1e-4 ? (tcltRad / tcltMax) : vec3(1.0);
			float tcltW = saturate(max(tcltConf, getLuminance(tcltRad) * 2.0)) * saturate(lightmap.x * 1.2);
			artificialColor = mix(artificialColor, tcltColor, tcltW * TCLT_TINT_STRENGTH);

			vec3 tcltGI = tcltRad * BRDF_D * TCLT_GI_STRENGTH;
			tcltGI *= localLightShadow(viewPos1.xyz, normalV, lightmap.x);
			L2 += tcltGI;

			tcltDbgRad = tcltRad;
			tcltDbgConf = max(tcltConf, tcltLocalConf);
			tcltDbgColor = tcltColor;
		#endif
		vec3 artificial = lightmap.x * artificialColor * (1. + GLOWING_BRIGHTNESS * glowingB) * diffuse;
		artificial += saturate(materialParams.emissiveness - lightmap.x) * diffuse * EMISSIVENESS_BRIGHTNESS;
		artificial *= localLightShadow(viewPos1.xyz, normalV, lightmap.x);
		
		if(lightningBolt > 0.5) color.rgb = vec3(1.0, 0.0, 0.0);
		// artificial += 1 * lightningBolt;

		

		
		
		
		#ifdef DISABLE_LEAKAGE_REPAIR
			lightMask = 1.0;
		#endif
		color.rgb = albedo * 0.01;
		color.rgb += skyLight * SKY_LIGHT_BRIGHTNESS;
		color.rgb += nightVision * diffuse * NIGHT_VISION_BRIGHTNESS;
		color.rgb += L2 * lightMask * RSM_BRIGHTNESS;
		color.rgb *= ao;
		color.rgb += direct * lightMask * DIRECT_LUMINANCE;

		if(isEyeInWater == 1){
			vec3 underWaterTransmit = saturate(exp(-(vec3(1.0) - waterFogColor) * (1.25 - lightmap.y) * 3));
			color.rgb *= underWaterTransmit * 1.0;
		}
		
		color.rgb += artificial;
		// color.rgb = sunColor * gi.rgb + lightmap.y * mix(sunColor, skyColor, SUN_SKY_BLEND - 0.05 * noRSM * lightmap.y) * mix(1.0, UoN * 0.5 + 0.5, 0.75);
		// color.rgb *= vec3(ao);

	}else{
		float d_p2a = RaySphereIntersection(earthPos, worldDir, vec3(0.0), earth_r + atmosphere_h).y;
		float d_p2e = RaySphereIntersection(earthPos, worldDir, vec3(0.0), earth_r).x;
		float d = d_p2e > 0.0 ? d_p2e : d_p2a;
		float dist1 = skyB > 0.5 ? d : worldDis1;

		float cloudTransmittance = 1.0;
		vec3 cloudScattering = vec3(0.0);
		float cloudHitLength = clamp(intersectHorizontalPlane(camera, worldDir, 650), 0.0, 20000.0);

		#ifdef VOLUMETRIC_CLOUDS
			float cloudScale = getCloudRenderScale();
			vec2 cloudOffset = vec2(1.0 - cloudScale, 0.0);
			vec2 cloud_uv = texcoord * cloudScale + cloudOffset;
			vec2 cloudMin = cloudOffset + invViewSize * 1.5;
			vec2 cloudMax = cloudOffset + vec2(cloudScale) - invViewSize * 1.5;
			if(!outScreen(cloud_uv, cloudMin, cloudMax) && camera.y < 5000.0)	{
				vec4 CT1_c = sampleCloudUpscaled(cloud_uv);
				if(dot(CT1_c.rgb, CT1_c.rgb) <= 1e-9){
					CT1_c.a = 1.0;
				}
				cloudScattering = CT1_c.rgb;
				cloudTransmittance = CT1_c.a;
			}
		#endif

		vec3 skyBaseColor = texture(colortex1, texcoord * 0.5 + 0.5).rgb * SKY_BASE_COLOR_BRIGHTNESS;
		vec3 celestial = drawCelestial(worldDir, cloudTransmittance, true);

		color.rgb = skyBaseColor;	
		color.rgb += celestial;
		cloudTransmittance = max(cloudTransmittance, 0.0);
		cloudScattering = max(cloudScattering, vec3(0.0));
		color.rgb = color.rgb * cloudTransmittance + cloudScattering;

		float VoL = saturate(dot(worldDir, sunWorldDir));
		float phase = saturate(hgPhase1(VoL, 0.66 - 0.56 * rainStrength));
		float crepuscularLight = 0.0;
		#ifdef CREPUSCULAR_LIGHT
			if(phase > 0.01 && sunRiseSetS + isNoonS > 0.001) crepuscularLight = computeCrepuscularLight(viewPos1) * phase;
		#endif
		if(cloudTransmittance < 1.0){

			color.rgb = 
				mix((skyBaseColor + celestial), color.rgb, 
					saturate(
						// mix(saturate(pow(getLuminance(cloudScattering), 1.0 - 0.45 * phase * sunRiseSetS)), 
						// 	exp(-cloudHitLength / (CLOUD_FADE_DISTANCE * (1.0 + 1.0 * phase * sunRiseSetS))) * 1.0, 
						// 	0.66)
						pow(exp(-cloudHitLength / (CLOUD_FADE_DISTANCE * (1.0 + 1.0 * phase * sunRiseSetS))), 
								remapSaturate(1.0 - saturate(getLuminance(cloudScattering - 0.5) + 0.05), 0.0, 1.0, 1.0, 2.0))
					)
				);
			
		}
		color.rgb += pow(crepuscularLight, 1.0) * sunColor * max3(0.6 * sunRiseSetS, 5.0 * rainStrength, 0.05 * isNoonS) * saturate(1.0 - isNightS)
					* remapSaturate(camera.y, 600.0, 1000.0, 1.0, 0.0);
		// color.rgb = vec3(computeCrepuscularLight(viewPos1));
		// color.rgb = vec3(crepuscularLight);
	}
	
	color.rgb = max(BLACK, color.rgb);
	#if DEBUG_TCLT == 1
		color.rgb = tcltDbgRad * 4.0 + vec3(tcltDbgConf);
	#elif DEBUG_TCLT == 2
		color.rgb = tcltDbgRad;
	#elif DEBUG_TCLT == 3
		color.rgb = tcltDbgColor;
	#elif DEBUG_TCLT == 4
		color.rgb = vec3(1.0, 0.0, 1.0);
	#endif
	// color.rgb = vec3(1.0 - texture(colortex3, texcoord * 0.5).a);
	
	// if(dhTerrain > 0.5) color.rgb = vec3(1.0 - texture(colortex1, texcoord * 0.5).a);
	
	CT4.rg = pack4x8To2x16(vec4(albedo, ao));

/* DRAWBUFFERS:0249 */
	gl_FragData[0] = color;
	gl_FragData[1] = CT2;
	gl_FragData[2] = CT4;
	gl_FragData[3] = vec4(velocity, 0.0, 1.0);
}

#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////BY ZY//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef VSH

void main() {
	sunViewDir = normalize(sunPosition);
	moonViewDir = normalize(moonPosition);
	lightViewDir = normalize(shadowLightPosition);

	sunWorldDir = normalize(viewPosToWorldPos(vec4(sunPosition, 0.0)).xyz);
    moonWorldDir = normalize(viewPosToWorldPos(vec4(moonPosition, 0.0)).xyz);
    lightWorldDir = normalize(viewPosToWorldPos(vec4(shadowLightPosition, 0.0)).xyz);

	isNoon = saturate(dot(sunWorldDir, upWorldDir) * NOON_DURATION);
	isNight = saturate(dot(moonWorldDir, upWorldDir) * NIGHT_DURATION);
	sunRiseSet = saturate(1 - isNoon - isNight);

	isNoonS = saturate(dot(sunWorldDir, upWorldDir) * NOON_DURATION_SLOW);
	isNightS = saturate(dot(moonWorldDir, upWorldDir) * NIGHT_DURATION_SLOW);
	sunRiseSetS = saturate(1 - isNoonS - isNightS);

	sunColor = getSunColor();
	skyColor = getSkyColor();

	gl_Position = ftransform();
	texcoord = (gl_TextureMatrix[0] * gl_MultiTexCoord0).xy;
}

#endif
