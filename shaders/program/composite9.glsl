varying vec2 texcoord;

#include "/lib/uniform.glsl"
#include "/lib/settings.glsl"
#include "/lib/common/utils.glsl"
#include "/lib/common/position.glsl"

#ifdef FSH

const bool shadowtex0Mipmap = false;
const bool shadowtex1Mipmap = false;
const bool shadowcolor0Mipmap = false;
const bool shadowcolor1Mipmap = false;

void main() {
    vec4 prev = texture(colortex10, texcoord);

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

    vec4 outColor = prev;
    if (best >= LOCAL_LIGHT_THRESHOLD) {
        vec4 viewPos = screenPosToViewPos(vec4(unTAAJitter(bestUV), bestDepth, 1.0));
        vec4 cur = vec4(viewPos.xyz, best);
        outColor = mix(prev, cur, LOCAL_LIGHT_HISTORY);
    } else {
        outColor = prev * (1.0 - LOCAL_LIGHT_DECAY);
    }

/* DRAWBUFFERS:A */
    gl_FragData[0] = outColor;
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
