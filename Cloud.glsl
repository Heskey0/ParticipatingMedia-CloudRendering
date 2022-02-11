
mat3 m = mat3( 0.00,  0.80,  0.60,
              -0.80,  0.36, -0.48,
              -0.60, -0.48,  0.64 );
float hash( float n )
{
    return fract(sin(n)*43758.5453);
}

//梯度噪声产生的纹理具有连续性，所以经常用来模拟山脉、云朵等具有连续性的物质，该类噪声的典型代表是Perlin Noise。
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);    //老版一维perlin noise  
    //https://blog.csdn.net/candycat1992/article/details/50346469
    //https://zhuanlan.zhihu.com/p/68507311

    float n = p.x + p.y*57.0 + 113.0*p.z;

    float res = mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
    return res;
}

//Fractal Brownian Motion
//分型布朗运动
float fbm( vec3 p )
{
    float f;
    f  = 0.5000*noise( p ); p = m*p*2.02;
    f += 0.2500*noise( p ); p = m*p*2.03;
    f += 0.12500*noise( p ); p = m*p*2.01;
    f += 0.06250*noise( p );
    return f;
}
/////////////////////////////////////

float stepUp(float t, float len, float smo)
{
  float tt = mod(t += smo, len);
  float stp = floor(t / len) - 1.0;
  return smoothstep(0.0, smo, tt) + stp;
}

//polynomial smooth minimum
//多项式平滑最小值
//https://zhuanlan.zhihu.com/p/281600585
// iq's smin
float smin( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}


float map( in vec3 p )
{
	vec3 q = p - vec3(0.0,0.5,1.0)*iTime;
    float f = fbm(q);

    float torus = 1. - sdTorus(p * 2.0, vec2(6.0, 0.005)) + f * 3.5;
    float s1 = 1.0 - length(p * vec3(0.5, 1.0, 0.5)) + f * 2.2;
    float s2 = 1.0 - length(p * vec3(0.1, 1.0, 0.2)) + f * 2.5;
    float s3 = 1.0 - smin(smin(
                           length(p * 1.0 - vec3(cos(iTime * 3.0) * 6.0, sin(iTime * 2.0) * 5.0, 0.0)),
                           length(p * 2.0 - vec3(0.0, sin(iTime) * 4.0, cos(iTime * 2.0) * 3.0)), 4.0),
                           length(p * 3.0 - vec3(cos(iTime * 2.0) * 3.0, 0.0, sin(iTime * 3.3) * 7.0)), 4.0) + f * 2.5;

    return min(max(0.0,torus), 1.0);
    

}

float jitter;

#define MAX_STEPS 48
#define SHADOW_STEPS 8
#define VOLUME_LENGTH 15.
#define SHADOW_LENGTH 2.

// https://shaderbits.com/blog/creating-volumetric-ray-marcher
vec4 cloudMarch(vec3 p, vec3 ray)
{

    float stepLength = VOLUME_LENGTH / float(MAX_STEPS);
    float shadowStepLength = SHADOW_LENGTH / float(SHADOW_STEPS);
    vec3 light = normalize(vec3(1.0, 2.0, 1.0));

    vec4 sum = vec4(0., 0., 0., 1.);
    
    vec3 pos = p + ray * jitter * stepLength;
    
    for (int i = 0; i < MAX_STEPS; i++)
    {
        if (sum.a < 0.1) {
        	break;
        }
        float d = map(pos);     //sdf
    
        if( d > 0.001)
        {
            vec3 lpos = pos + light * jitter * shadowStepLength;
            float shadow = 0.;
    
            for (int s = 0; s < SHADOW_STEPS; s++)
            {
                lpos += light * shadowStepLength;
                float lsample = map(lpos);
                shadow += lsample;
            }
    
            //totalDensity
            // = sum(sampleDensity(lpos) * stepSize)
            float density = clamp(d * VOLUME_LENGTH / float(MAX_STEPS), 0.0, 1.0);
            //density = map()
            density = sin(iTime * .3 - 1.7) * 1. + 0.9;

            //transmittance
            // = exp(density1 * lightAbsorption)
            //  *exp(density2 * lightAbsorption)
            //  *...
            float s = exp(-shadow * SHADOW_LENGTH / float(SHADOW_STEPS) * 1.0);

            sum.rgb += vec3(s * density) * vec3(1.1, 0.9, .5) * sum.a;
            //sum.rgb *= s;
            sum.a *= 1.-density;

            sum.rgb += exp(-map(pos + vec3(0,0.25,0.0)) * .2) * density * vec3(0.15, 0.45, 1.1) * sum.a;
        }
        pos += ray * stepLength;
    }
 
    //return sum * vec4(1.1, 0.9, .5, 1);
    return sum;
}

mat3 camera(vec3 ro, vec3 ta, float cr )
{
	vec3 cw = normalize(ta - ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 p = (fragCoord.xy * 2.0 - iResolution.xy) / min(iResolution.x, iResolution.y); //<-1, 1>
    //TODO:replace one line
    //jitter = hash(p.x + p.y * 57.0 + iTime);
    jitter = 0.1;
    
    vec3 ro = vec3(cos(iTime * .333) * 8.0, -5.5, sin(iTime * .333) * 8.0);
    vec3 ta = vec3(0.0, 1., 0.0);
    mat3 c = camera(ro, ta, 0.0);
    vec3 ray = c * normalize(vec3(p, 1.75));
    vec4 col = cloudMarch(ro, ray);
    vec3 result = col.rgb + mix(vec3(0.3, 0.6, 1.0), vec3(0.05, 0.35, 1.0), p.y + 0.75) * (col.a);
    
    float sundot = clamp(dot(ray,normalize(vec3(1.0, 2.0, 1.0))),0.0,1.0);
    result += 0.4*vec3(1.0,0.7,0.3)*pow( sundot, 4.0 );

    result = pow(result, vec3(1.0/2.2));
    //result = col.rgb;
    
    fragColor = vec4(result,1.0);
}