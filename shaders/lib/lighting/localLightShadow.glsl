#ifndef LOCAL_LIGHT_SHADOW_GLSL
#define LOCAL_LIGHT_SHADOW_GLSL

float localLightShadow(vec3 viewPos, vec3 normalV, float blockLight) {
    #ifdef LOCAL_LIGHT_SHADOWS
        if (blockLight < 0.05) return 1.0;
        vec4 lp = texture(colortex10, texcoord);
        if (lp.a < LOCAL_LIGHT_THRESHOLD) return 1.0;

        vec3 toLight = lp.rgb - viewPos;
        float dist = length(toLight);
        if (dist < 0.001 || dist > LOCAL_LIGHT_MAX_DIST) return 1.0;

        vec3 dir = toLight / dist;
        float ndotl = saturate(dot(normalV, dir));
        if (ndotl <= 0.0) return 1.0;

        float shadow = 1.0;
        float stepLen = dist / float(LOCAL_LIGHT_SHADOW_STEPS);
        for (int i = 1; i < LOCAL_LIGHT_SHADOW_STEPS; ++i) {
            vec3 p = viewPos + dir * (stepLen * float(i));
            vec3 sp = viewPosToScreenPos(vec4(p, 1.0)).xyz;
            if (outScreen(sp.xy)) {
                shadow = 1.0;
                break;
            }
            float depth = texture(depthtex1, sp.xy).r;
            if (depth < sp.z - LOCAL_LIGHT_SHADOW_BIAS) {
                shadow = 0.0;
                break;
            }
        }
        return mix(1.0, shadow, ndotl);
    #else
        return 1.0;
    #endif
}

#endif
