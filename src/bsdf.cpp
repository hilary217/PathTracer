#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface

    *pdf = 1.0;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  double tan = sin_theta(h) / cos_theta(h);
  return exp(-pow(tan, 2) / pow(alpha, 2)) / (PI * pow(alpha, 2) * pow(cos_theta(h), 4));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.

  Spectrum Rs, Rp;
  double c = cos_theta(wi);
  double csquare = pow(c, 2);
  Spectrum s2 = (eta * eta) + (k * k);
  Rs = (s2 - 2 * eta * c + csquare) / (s2 + 2 * eta * c + csquare);
  Rp = (s2 * csquare - 2 * eta * c + 1) / (s2 * csquare + 2 * eta * c + 1);
  return (Rs + Rp) / 2;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here
  Vector3D n(0,0,1);
  if((dot(n, wo) > 0) && (dot(n, wi) > 0)){
      Vector3D h = wo + wi;
      h.normalize();
      return F(wi) * G(wo, wi) * D(h) / (4.0 * dot(n, wo) * dot(n, wi));
  }else{
      return Spectrum();
  }
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

    //*wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
    //return MicrofacetBSDF::f(wo, *wi);

    Vector2D xy = sampler.get_sample();


    float theta= atan(sqrt(-pow(alpha, 2) * log(1.0 - xy.x)));
    float phi = 2.0 * PI * xy.y;
    Vector3D h(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    *wi = -wo + 2.0 * dot(wo, h) * h;

    float ptheta = 2.0 * sin(theta) * exp(- pow(tan(theta), 2) / pow(alpha, 2)) / (pow(alpha, 2) * pow(cos(theta), 3));
    float pphi = 1.0 / (2.0 * PI);
    float pwh = ptheta * pphi / sin(theta);
    float pwi = pwh / (4.0 * dot(*wi, h));
    *pdf = pwi;
    return MicrofacetBSDF::f(wo, *wi);


}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.
    if(!refract(wo,wi,ior)){
        *pdf = 1.0;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    else{
        // Compute Schlick's reflection coefficient R
        float n1 = 1.0;
        float n2 = ior;
        float R_o = pow(((n1-n2)/ (n1+n2)), 2);
        float R = R_o + (1.0 - R_o) * pow((1.0 - abs_cos_theta(wo)), 5);
        if(coin_flip(R)){
            *pdf = R;
            reflect(wo, wi);
            return R * reflectance / abs_cos_theta(*wi);
        }
        else{
            refract(wo,wi,ior);
            *pdf = 1.0 - R;
            float eta;
            if(wo.z >= 0)
                eta = 1.0 / ior;
            else
                eta = ior;

            return (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
        }
    }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);


}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  float eta;
  float wi_x_p, wi_y_p, wi_z_p;
  if(wo.z >= 0){
     eta = 1.0/ior;
     float costheta = wo.z;
     float is = 1 - pow(eta,2)*(1-pow(costheta,2));
     if(is<0){
         return false;
     }else{
         wi_z_p = - sqrt(is);
         wi_x_p = - eta* wo.x;
         wi_y_p = - eta * wo.y;
     }

  }
  else{
      eta = ior;
      float costheta = wo.z;
      float is = 1 - pow(eta,2)*(1-pow(costheta,2));
      if(is<0){
          return false;
      }else{
          wi_z_p =  sqrt(is);
          wi_x_p = - eta* wo.x;
          wi_y_p = - eta * wo.y;
      }

  }



  *wi = Vector3D(wi_x_p, wi_y_p, wi_z_p);

  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
