varying vec2 texcoord;

#include "/lib/uniform.glsl"
#include "/lib/settings.glsl"
#include "/lib/common/utils.glsl"
#include "/lib/common/position.glsl"
#include "/lib/common/normal.glsl"
#include "/lib/common/noise.glsl"
#include "/lib/camera/colorToolkit.glsl"

#ifdef FSH

const bool shadowtex0Mipmap = false;
const bool shadowtex1Mipmap = false;
const bool shadowcolor0Mipmap = false;
const bool shadowcolor1Mipmap = false;

void main() {
    vec4 prevLight = texture(colortex10, texcoord);

    float best = 0.0;
    float bestDepth = 1.0;
    vec2 bestUV = texcoord;
    vec2 step = LOCAL_LIGHT_TILE_SIZE * invViewSize * 0.5;

    for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
            vec2 uv = texcoord + vec2(float(i), float(j)) * step;
            if (outScreen(uv)) continue;

            float depth = texture(depthtex1, uv).r;
            if (depth >= 1.0) continue;

            vec2 lm = texture(colortex5, uv).ba;
            float blockL = lm.x;
            if (blockL > best) {
                best = blockL;
                bestDepth = depth;
                bestUV = uv;
            }
        }
    }

    vec4 localOut = prevLight;
    if (best >= LOCAL_LIGHT_THRESHOLD) {
        vec4 viewPos = screenPosToViewPos(vec4(unTAAJitter(bestUV), bestDepth, 1.0));
        vec4 cur = vec4(viewPos.xyz, best);
        localOut = mix(prevLight, cur, LOCAL_LIGHT_HISTORY);
    } else {
        localOut = prevLight * (1.0 - LOCAL_LIGHT_DECAY);
    }

    vec4 tcltOut = texture(colortex11, texcoord);
    #ifdef TCLT_COLORED_LIGHTING
        vec2 cacheSize = vec2(textureSize(colortex11, 0));
        vec2 tileUv = (gl_FragCoord.xy + 0.5) / cacheSize;
        vec2 tileStep = 1.0 / cacheSize;

        vec2 jitter = hash42(gl_FragCoord.xy).rg * 2.0 - 1.0;
        jitter *= 0.35;

        vec3 sum = vec3(0.0);

        vec2 offsets[4];
        offsets[0] = vec2( 0.25,  0.25);
        offsets[1] = vec2(-0.25,  0.25);
        offsets[2] = vec2( 0.25, -0.25);
        offsets[3] = vec2(-0.25, -0.25);

        for (int s = 0; s < 4; ++s) {
            vec2 uv = tileUv + (offsets[s] + jitter) * tileStep;
            if (outScreen(uv)) continue;

            float depth = texture(depthtex1, uv).r;
            if (depth >= 1.0) continue;

            vec3 albedo = pow(texture(colortex0, uv).rgb, vec3(2.2));
            vec2 lm = texture(colortex5, uv).ba;
            float blockL = lm.x;

            vec4 spec = unpack2x16To4x8(texture(colortex4, uv).ba);
            float emissive = spec.a;
            if (emissive * 255.0 > 254.1) emissive = 0.0;

            float energy = max(emissive * TCLT_EMISSIVE_GAIN, blockL * TCLT_BLOCKLIGHT_GAIN);
            energy = max(0.0, energy - TCLT_EMISSIVE_TH);
            float sat = max3(albedo.r, albedo.g, albedo.b);
            energy *= saturate(sat * 1.2);

            sum += albedo * energy;
        }

        vec3 cur = sum * 0.25;

        float depthC = texture(depthtex1, tileUv).r;
        vec3 normalC = normalize(normalDecode(texture(colortex5, tileUv).rg));
        float depthCLin = linearizeDepth(depthC);

        vec3 neighSum = vec3(0.0);
        float neighW = 0.0;

        vec2 nb[4];
        nb[0] = vec2(1.0, 0.0);
        nb[1] = vec2(-1.0, 0.0);
        nb[2] = vec2(0.0, 1.0);
        nb[3] = vec2(0.0, -1.0);

        for (int i = 0; i < 4; ++i) {
            vec2 nuv = tileUv + nb[i] * tileStep;
            if (outScreen(nuv)) continue;

            float nd = texture(depthtex1, nuv).r;
            if (nd >= 1.0) continue;

            float ndLin = linearizeDepth(nd);
            float depthW = step(abs(ndLin - depthCLin), TCLT_DEPTH_TH);
            vec3 nn = normalize(normalDecode(texture(colortex5, nuv).rg));
            float normalW = smoothstep(TCLT_NORMAL_TH, 1.0, saturate(dot(nn, normalC)));
            float w = depthW * normalW;
            if (w <= 0.0) continue;

            vec3 nc = texture(colortex11, nuv).rgb;
            neighSum += nc * w;
            neighW += w;
        }

        vec3 neigh = neighW > 1e-4 ? (neighSum / neighW) : vec3(0.0);
        vec3 curProp = mix(cur, neigh, TCLT_PROPAGATION);

        vec3 prev = texture(colortex11, tileUv).rgb;
        float valid = 1.0;
        vec2 prevUv = tileUv;
        if (depthC >= 1.0) valid = 0.0;
        if (valid > 0.0) {
            vec4 viewPos = screenPosToViewPos(vec4(unTAAJitter(tileUv), depthC, 1.0));
            vec4 worldPos = viewPosToWorldPos(viewPos);
            prevUv = getPrePos(worldPos).xy;
            if (outScreen(prevUv)) valid = 0.0;
        }
        if (valid > 0.0) {
            float preDepth = texture(depthtex1, prevUv).r;
            float preLin = linearizeDepth(preDepth);
            float depthDiff = abs(preLin - depthCLin);
            valid = step(depthDiff, TCLT_DEPTH_TH);
            prev = texture(colortex11, prevUv).rgb;
        }

        vec3 blended = mix(prev, curProp, TCLT_HISTORY);
        float signal = saturate(getLuminance(curProp) * 4.0);
        vec3 outColor = mix(prev * (1.0 - TCLT_DECAY), blended, signal);
        float conf = saturate(getLuminance(outColor) * 2.0);
        tcltOut = vec4(outColor, conf);
    #endif
    #if DEBUG_TCLT == 1
        tcltOut = vec4(1.0, 0.0, 1.0, 1.0);
    #endif

/* DRAWBUFFERS:AB */
    gl_FragData[0] = localOut;
    gl_FragData[1] = tcltOut;
}

#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////BY ZYPanDa/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef VSH

void main() {
    gl_Position = ftransform();
    texcoord = (gl_TextureMatrix[0] * gl_MultiTexCoord0).xy;
}

#endif
