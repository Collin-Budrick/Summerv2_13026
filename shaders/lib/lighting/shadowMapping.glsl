float simpleShadowMapping(vec4 worldPos){
    vec4 shadowPos = getShadowPos(worldPos);
    float shadow = texture(shadowtex0, vec3(shadowPos.xy, shadowPos.z - 0.0005)).r;

    return shadow;
}

float PenumbraMask(vec2 uv){
    if(outScreen(uv)) return 0.0;

    vec2 offset[4] = {
        vec2(-1.0, 1.0),
        vec2(1.0, 1.0),
        vec2(-1.0, -1.0),
        vec2(1.0, -1.0)
    };

    float shadow = 0.0;
    for(int i = 0; i < 4; i++){
        vec2 curUV = uv + offset[i] * invViewSize;
        float depth = texture(depthtex1, curUV).r;
        vec4 viewPos = screenPosToViewPos(vec4(unTAAJitter(curUV), depth, 1.0));
	    vec4 worldPos = viewPosToWorldPos(viewPos);

        float curShadow = simpleShadowMapping(worldPos);
        shadow += curShadow;
    }
    shadow *= 0.25;

    return (shadow < 1.0 && shadow > 0.0) ? 1.0 : 0.0;
    // return shadow;
}

float PenumbraMaskBlur(vec2 uv, vec2 dir){
    if(outScreen(uv * 4.0)) return 0.0;
    float depth1 = texture(depthtex1, uv * 4.0).r;
    vec4 viewPos1 = screenPosToViewPos(vec4(unTAAJitter(uv * 4.0), depth1, 1.0));
	vec4 worldPos1 = viewPosToWorldPos(viewPos1);
	float worldDis1 = length(worldPos1);

    float r = 2.5;
    float scale = clamp(1.0 / (worldDis1 * 0.08), 1.0, 20.0);
    float PenumbraMaskBlur = 0.0;
    for(float i = -r; i < r + 0.1; i++){
        vec2 curUV = uv + i * dir * scale;
        PenumbraMaskBlur = max(PenumbraMaskBlur, textureLod(colortex1, curUV, log(scale)).r);
    }
    return step(0.001, PenumbraMaskBlur);
}


#ifdef FSH
float blockerSearch(sampler2D shadowMap, vec3 shadowPos, float radius, float quality){
    float c = 0.0;
    float blocker = 0.0;

    float radiusStep = radius;

    float noise = temporalBayer64(gl_FragCoord.xy);
    float firstAngle = noise * _2PI * 17.3333333;
    vec2 curDir = vec2(cos(firstAngle), sin(firstAngle));

    float rotAngle = GOLDEN_ANGLE;
    float sinRot = sin(rotAngle);
    float cosRot = cos(rotAngle);
    mat2 rotMatrix = mat2(cosRot, -sinRot, sinRot, cosRot);

    radius *= noise;

    for (int i = 0; i < quality; i++) {
        vec2 offset = curDir * pow(radius, 0.75);

        float dBlocker = textureLod(shadowMap, shadowPos.xy + offset / shadowMapResolution, 0.0).r;
        if(shadowPos.z > dBlocker){
            blocker += dBlocker;
            c++;
        }

        radius += radiusStep;
        curDir = rotMatrix * curDir;
    }
    if(c <= 0.8) return -1.0;

    return blocker / c;
}

float PCF(sampler2DShadow shadowMap, vec3 shadowPos, float radius, float quality){
    float c = 0.0;
    float shade = 0.0;

    float radiusStep = radius;

    float noise = temporalBayer64(gl_FragCoord.xy);
    float firstAngle = noise * _2PI * 17.3333333;
    vec2 curDir = vec2(cos(firstAngle), sin(firstAngle));

    float rotAngle = GOLDEN_ANGLE;
    float sinRot = sin(rotAngle);
    float cosRot = cos(rotAngle);
    mat2 rotMatrix = mat2(cosRot, -sinRot, sinRot, cosRot);

    radius *= noise;

    for (int i = 0; i < quality; i++) {
        vec2 offset = curDir * pow(radius, 0.75);

        shade += textureLod(shadowMap, vec3(shadowPos.xy + offset / shadowMapResolution, shadowPos.z), 0.0).r;
        c++;

        radius += radiusStep;
        curDir = rotMatrix * curDir;
    }

    return shade / c;
}
#endif

#ifndef GBF
#ifdef SHADOW_RESERVOIR
float stableBlueNoise(vec2 fragCoord){
    return texelFetch(noisetex, ivec2(fragCoord) % noiseTextureResolution, 0).r;
}

float blockerSearchStable(sampler2D shadowMap, vec3 shadowPos, float radius, float quality){
    float c = 0.0;
    float blocker = 0.0;

    float radiusStep = radius;

    float noise = stableBlueNoise(gl_FragCoord.xy);
    float firstAngle = noise * _2PI * 17.3333333;
    vec2 curDir = vec2(cos(firstAngle), sin(firstAngle));

    float rotAngle = GOLDEN_ANGLE;
    float sinRot = sin(rotAngle);
    float cosRot = cos(rotAngle);
    mat2 rotMatrix = mat2(cosRot, -sinRot, sinRot, cosRot);

    radius *= noise;

    for (int i = 0; i < int(quality); i++) {
        vec2 offset = curDir * pow(radius, 0.75);

        float dBlocker = textureLod(shadowMap, shadowPos.xy + offset / shadowMapResolution, 0.0).r;
        if(shadowPos.z > dBlocker){
            blocker += dBlocker;
            c++;
        }

        radius += radiusStep;
        curDir = rotMatrix * curDir;
    }
    if(c <= 0.8) return -1.0;

    return blocker / c;
}

bool validateShadowReservoir(vec2 uv, float curDepthLin, vec3 curNormalV){
    if(outScreen(uv)) return false;
    float depth = texture(depthtex1, uv).r;
    float depthLin = linearizeDepth(depth);
    if(abs(depthLin - curDepthLin) > SHADOW_RESERVOIR_DEPTH_TH * max(curDepthLin, 1.0)) return false;

    vec2 nEnc = texelFetch(colortex5, ivec2(uv * viewSize), 0).rg;
    vec3 nV = normalize(normalDecode(nEnc));
    return dot(nV, curNormalV) > SHADOW_RESERVOIR_NORMAL_TH;
}

float estimatePenumbra(vec3 shadowPos, float sssWrap){
    float dReceiver = shadowPos.z;
    float dBlocker = blockerSearchStable(shadowtex1, shadowPos.xyz, 2.5 * shadowMapScale, SHADOW_RESERVOIR_BLOCKER_SAMPLES);
    if(dBlocker <= 0.0) return 0.0;

    float penumbra = max(0.0, (dReceiver - dBlocker) / max(dBlocker, 1e-4));
    penumbra = penumbra * SHADOW_ANGULAR_RADIUS * 100.0;
    penumbra += 0.5 * sssWrap;
    penumbra = min(penumbra, SHADOW_PENUMBRA_MAX);
    return penumbra;
}

float shadowMapping(vec4 worldPos, vec3 normal, float sssWrap, float packedIn, out float packedOut){
    #ifndef DISTANT_HORIZONS
        if(skyB > 0.5) return 1.0;
    #endif

    float worldDis = length(worldPos.xyz);
    vec4 shadowPos = getShadowPos(worldPos);

    float disFactor = saturate(worldDis / shadowDistance);
    float dirFactor = abs(dot(normal, lightWorldDir));

    float offset = 0.05;
    if(plants < 0.5) {
        offset = 0.05 + 1.0 * disFactor * (1 - dirFactor);
    }
    worldPos.xyz += normal * offset;
    shadowPos = getShadowPos(worldPos);

    float penumbra = estimatePenumbra(shadowPos.xyz, sssWrap);
    if(plants > 0.5){
        penumbra = 1.5;
    }
    shadowPos.z -= 0.00005;
    penumbra *= saturate(saturate(1.0 - disFactor) + 0.1);

    float baseRadius = (0.25 + penumbra + rainStrength) * shadowMapScale * SHADOW_SOFTNESS;
    float centerVis = textureLod(shadowtex0, vec3(shadowPos.xy, shadowPos.z), 0.0).r;

    float sumW = 0.0;
    float chosenAngle = 0.0;

    float noise = stableBlueNoise(gl_FragCoord.xy);
    float angle = noise * _2PI;
    float rotAngle = GOLDEN_ANGLE;

    for(int i = 0; i < SHADOW_RESERVOIR_CANDIDATES; ++i){
        vec2 dir = vec2(cos(angle), sin(angle));
        float jitter = fract(noise + float(i) * 0.6180339);
        float radius = baseRadius * (0.25 + 0.75 * jitter);
        vec2 offset = dir * (radius / shadowMapResolution);
        float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
        float w = max(abs(vis - centerVis), 0.1);
        sumW += w;

        float pick = fract(noise + float(i) * 0.37);
        if(pick * sumW < w){
            chosenAngle = angle;
        }

        angle += rotAngle;
    }

    vec3 shadowReflData = unpackShadowReflData(packedIn);
    float reflN = shadowReflData.z;
    float curDepthLin = linearizeDepth(texture(depthtex1, texcoord).r);
    vec3 curNormalV = normalize(mat3(gbufferModelView) * normal);

    vec3 prePos = getPrePos(worldPos);
    if(validateShadowReservoir(prePos.xy, curDepthLin, curNormalV)){
        vec3 preShadowRefl = unpackShadowReflData(texture(colortex2, prePos.xy).a);
        float prevAngle = preShadowRefl.x * _2PI;
        float prevWeight = preShadowRefl.y * SHADOW_RESERVOIR_MAX * SHADOW_RESERVOIR_HISTORY;
        if(prevWeight > 0.0){
            vec2 dir = vec2(cos(prevAngle), sin(prevAngle));
            vec2 offset = dir * (baseRadius / shadowMapResolution);
            float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
            float w = prevWeight;
            sumW += w;
            float pick = fract(noise + 0.73);
            if(pick * sumW < w){
                chosenAngle = prevAngle;
            }
        }
    }

    #ifdef SHADOW_RESERVOIR_SPATIAL
        ivec2 pix = ivec2(gl_FragCoord.xy);
        ivec2 p = pix + ivec2(-1, 0);
        vec2 uv = (vec2(p) + 0.5) * invViewSize;
        if(validateShadowReservoir(uv, curDepthLin, curNormalV)){
            vec3 nbShadowRefl = unpackShadowReflData(texelFetch(colortex2, p, 0).a);
            float nbAngle = nbShadowRefl.x * _2PI;
            float nbWeight = nbShadowRefl.y * SHADOW_RESERVOIR_MAX * 0.5;
            if(nbWeight > 0.0){
                vec2 dir = vec2(cos(nbAngle), sin(nbAngle));
                vec2 offset = dir * (baseRadius / shadowMapResolution);
                float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
                float w = nbWeight;
                sumW += w;
                float pick = fract(noise + 0.17);
                if(pick * sumW < w){
                    chosenAngle = nbAngle;
                }
            }
        }

        p = pix + ivec2(1, 0);
        uv = (vec2(p) + 0.5) * invViewSize;
        if(validateShadowReservoir(uv, curDepthLin, curNormalV)){
            vec3 nbShadowRefl = unpackShadowReflData(texelFetch(colortex2, p, 0).a);
            float nbAngle = nbShadowRefl.x * _2PI;
            float nbWeight = nbShadowRefl.y * SHADOW_RESERVOIR_MAX * 0.5;
            if(nbWeight > 0.0){
                vec2 dir = vec2(cos(nbAngle), sin(nbAngle));
                vec2 offset = dir * (baseRadius / shadowMapResolution);
                float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
                float w = nbWeight;
                sumW += w;
                float pick = fract(noise + 0.46);
                if(pick * sumW < w){
                    chosenAngle = nbAngle;
                }
            }
        }

        p = pix + ivec2(0, -1);
        uv = (vec2(p) + 0.5) * invViewSize;
        if(validateShadowReservoir(uv, curDepthLin, curNormalV)){
            vec3 nbShadowRefl = unpackShadowReflData(texelFetch(colortex2, p, 0).a);
            float nbAngle = nbShadowRefl.x * _2PI;
            float nbWeight = nbShadowRefl.y * SHADOW_RESERVOIR_MAX * 0.5;
            if(nbWeight > 0.0){
                vec2 dir = vec2(cos(nbAngle), sin(nbAngle));
                vec2 offset = dir * (baseRadius / shadowMapResolution);
                float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
                float w = nbWeight;
                sumW += w;
                float pick = fract(noise + 0.71);
                if(pick * sumW < w){
                    chosenAngle = nbAngle;
                }
            }
        }

        p = pix + ivec2(0, 1);
        uv = (vec2(p) + 0.5) * invViewSize;
        if(validateShadowReservoir(uv, curDepthLin, curNormalV)){
            vec3 nbShadowRefl = unpackShadowReflData(texelFetch(colortex2, p, 0).a);
            float nbAngle = nbShadowRefl.x * _2PI;
            float nbWeight = nbShadowRefl.y * SHADOW_RESERVOIR_MAX * 0.5;
            if(nbWeight > 0.0){
                vec2 dir = vec2(cos(nbAngle), sin(nbAngle));
                vec2 offset = dir * (baseRadius / shadowMapResolution);
                float vis = textureLod(shadowtex0, vec3(shadowPos.xy + offset, shadowPos.z), 0.0).r;
                float w = nbWeight;
                sumW += w;
                float pick = fract(noise + 0.93);
                if(pick * sumW < w){
                    chosenAngle = nbAngle;
                }
            }
        }
    #endif

    float totalW = min(sumW, SHADOW_RESERVOIR_MAX);
    float angleN = sumW > 0.0 ? fract(chosenAngle / _2PI) : 0.0;
    float shadowWN = totalW / SHADOW_RESERVOIR_MAX;
    packedOut = packShadowReflData(vec3(angleN, shadowWN, reflN));

    vec2 baseDir = vec2(cos(chosenAngle), sin(chosenAngle));
    vec2 baseOffset = baseDir * (baseRadius / shadowMapResolution);
    float shade = 0.0;
    float count = 0.0;
    float tapRadius = max(0.75, baseRadius * 0.5);
    vec2 offsetX = vec2(tapRadius / shadowMapResolution, 0.0);
    vec2 offsetY = vec2(0.0, tapRadius / shadowMapResolution);

    vec2 uv0 = shadowPos.xy + baseOffset;
    shade += textureLod(shadowtex0, vec3(uv0, shadowPos.z), 0.0).r;
    count += 1.0;

    vec2 uv1 = shadowPos.xy + baseOffset + offsetX;
    vec2 uv2 = shadowPos.xy + baseOffset - offsetX;
    shade += textureLod(shadowtex0, vec3(uv1, shadowPos.z), 0.0).r;
    shade += textureLod(shadowtex0, vec3(uv2, shadowPos.z), 0.0).r;
    count += 2.0;

    vec2 uv3 = shadowPos.xy + baseOffset + offsetY;
    vec2 uv4 = shadowPos.xy + baseOffset - offsetY;
    shade += textureLod(shadowtex0, vec3(uv3, shadowPos.z), 0.0).r;
    shade += textureLod(shadowtex0, vec3(uv4, shadowPos.z), 0.0).r;
    count += 2.0;

    vec2 uv5 = shadowPos.xy + baseOffset + (offsetX + offsetY);
    vec2 uv6 = shadowPos.xy + baseOffset + (offsetX - offsetY);
    vec2 uv7 = shadowPos.xy + baseOffset + (-offsetX + offsetY);
    vec2 uv8 = shadowPos.xy + baseOffset + (-offsetX - offsetY);
    shade += textureLod(shadowtex0, vec3(uv5, shadowPos.z), 0.0).r;
    shade += textureLod(shadowtex0, vec3(uv6, shadowPos.z), 0.0).r;
    shade += textureLod(shadowtex0, vec3(uv7, shadowPos.z), 0.0).r;
    shade += textureLod(shadowtex0, vec3(uv8, shadowPos.z), 0.0).r;
    count += 4.0;

    shade /= max(count, 1.0);

    return saturate(shade);
}
#else
float shadowMapping(vec4 worldPos, vec3 normal, float sssWrap){
    #ifndef DISTANT_HORIZONS
        if(skyB > 0.5) return 1.0;
    #endif

    float worldDis = length(worldPos.xyz);
    vec4 shadowPos = getShadowPos(worldPos);

    float dReceiver = shadowPos.z;
    float dBlocker = blockerSearch(shadowtex1, shadowPos.xyz, 4 * shadowMapScale, BLOCKER_SEARCH_SAMPLES);
    float penumbra = max(0.0, 100 * (dReceiver - dBlocker) / dBlocker);

    float disFactor = saturate(worldDis / shadowDistance);
    float dirFactor = abs(dot(normal, lightWorldDir));

    float offset = 0.05;
    if(plants < 0.5) {
        offset = 0.05 + 1.0 * disFactor * (1 - dirFactor);
    }
    worldPos.xyz += normal * offset;
    shadowPos = getShadowPos(worldPos);

    penumbra += 0.5 * sssWrap;
    if(plants > 0.5){
        penumbra = 1.5;
    }
    shadowPos.z -= 0.00005;
    penumbra *= saturate(saturate(1.0 - disFactor) + 0.1);

    float shade = PCF(shadowtex0, shadowPos.xyz, (0.25 + penumbra + rainStrength) * shadowMapScale * SHADOW_SOFTNESS, SHADOW_SAMPLES);

    return saturate(shade);
}
#endif

vec3 getColorShadow(vec3 shadowPos, float shadow){  
    const float N_SAMPLE = COLOR_SHADOW_SAMPLES;

    vec3 colorShadow = vec3(0.0);

    if(shadow < 0.9){
        // 每次循环跨越的距离
        float radiusStep = 0.5;
        // 设置起点（初始角度，距离），应用noise
        float noise = temporalBayer64(gl_FragCoord.xy);
        float firstAngle = noise * _2PI * 13.333333;
        float radius = radiusStep * noise;

        for (int i = 0; i < N_SAMPLE; i++) {
            float angle = firstAngle + rotationArray[i];
            vec2 offset = vec2(cos(angle), sin(angle)) * pow(radius, 0.75);
            vec2 uv = shadowPos.xy + offset / shadowMapResolution;

            vec4 SC1 = textureLod(shadowcolor1, uv, 0);
            SC1.rgb = SC1.rgb * 2.0 - 1.0;
            bool isTranslucent = length(SC1.rgb) < 0.1;
            if(!isTranslucent) {
                radius += radiusStep;
                continue;
            }

            float z_sample = textureLod(shadowtex1, uv, 0).r;
            if(shadowPos.z - 0.0001 > z_sample){
                radius += radiusStep;
                continue;
            }

            vec4 SC0 = textureLod(shadowcolor0, uv, 0);
            colorShadow += toLinearR(SC0.rgb);

            radius += radiusStep;
        }
    }
    return colorShadow / N_SAMPLE;
}
#endif




float Chebychev(float t, float mean, float variance){
    return variance / (variance + (t - mean) * (t - mean));
}

float variance(vec4 shadowPos, out float mean, out float meanX2){
    mean = textureLod(shadowcolor1, shadowPos.xy, 3).b;
    float mean2 = mean * mean;

    meanX2 = textureLod(shadowcolor1, shadowPos.xy, 3).a;

    return meanX2 - mean2;
}

float VSSM(vec4 worldPos){
    vec4 shadowPos = getShadowPos(worldPos);
    float mean = textureLod(shadowcolor1, shadowPos.xy, 3).b;
    float meanX2 = textureLod(shadowcolor1, shadowPos.xy, 3).a;
    float variance = meanX2 - mean * mean;
    
    float t = shadowPos.z;
    float N1 = Chebychev(t, mean, variance);
    float N2 = 1.0 - N1;
    float z_occ = (mean - N1 * t) / N2;

    float penumbra = 120 * (t - z_occ) / z_occ;

    float lod = log2(penumbra);
    mean = textureLod(shadowcolor1, shadowPos.xy, lod).b;
    meanX2 = textureLod(shadowcolor1, shadowPos.xy, lod).a;
    variance = meanX2 - mean * mean;

    return saturate(Chebychev(t, mean, variance));
}
