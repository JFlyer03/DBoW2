/**
 * File: ScoringObject.cpp
 * Date: November 2011
 * Author: Dorian Galvez-Lopez
 * Description: functions to compute bow scores 
 * License: see the LICENSE.txt file
 *
 */

#include <cfloat>
#include "TemplatedVocabulary.h"
#include "BowVector.h"

using namespace DBoW2;

// If you change the type of WordValue, make sure you change also the
// epsilon value (this is needed by the KL method)
const double GeneralScoring::LOG_EPS = log(DBL_EPSILON); // FLT_EPSILON

// Double-pointer traversal reduces the lookup overhead
double L1Scoring::score(const BowVector &v1, const BowVector &v2) const
{
  auto v1_it = v1.begin(), v2_it = v2.begin();
  const auto v1_end = v1.end(), v2_end = v2.end();
  double score = 0;

  while(v1_it != v1_end && v2_it != v2_end)
  {
    if(v1_it->first == v2_it->first)
    {
      const WordValue &vi = v1_it->second;
      const WordValue &wi = v2_it->second;
      score += fabs(vi - wi) - fabs(vi) - fabs(wi);
      ++v1_it; 
      ++v2_it;
    }
    else if(v1_it->first < v2_it->first) 
      ++v1_it;
    else 
      ++v2_it;
  }

  score = -score / 2.0;
  return score;
}

// ...existing code...
#include <execution> // C++17 parallel algorithms

double L2Scoring::score(const BowVector &v1, const BowVector &v2) const
{
  // parallel traversal to speed up scoring
  std::vector<double> partialScores;
  partialScores.reserve(std::max(v1.size(), v2.size()));
  
  std::for_each(std::execution::par, v1.begin(), v1.end(), [&](const auto &p1) {
    auto it = v2.find(p1.first);
    if(it != v2.end())
    {
      // skips pure zero-valued elements
      if(p1.second != 0.0 || it->second != 0.0)
      {
        partialScores.push_back(std::fabs(p1.second - it->second)
                                - std::fabs(p1.second)
                                - std::fabs(it->second));
      }
    }
  });

  double score = 0;
  for(const auto &val : partialScores)
  {
    score += val;
  }

  return -score / 2.0;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double ChiSquareScoring::score(const BowVector &v1, const BowVector &v2) 
  const
{
  BowVector::const_iterator v1_it, v2_it;
  const BowVector::const_iterator v1_end = v1.end();
  const BowVector::const_iterator v2_end = v2.end();
  
  v1_it = v1.begin();
  v2_it = v2.begin();
  
  double score = 0;
  
  // all the items are taken into account
  
  while(v1_it != v1_end && v2_it != v2_end)
  {
    const WordValue& vi = v1_it->second;
    const WordValue& wi = v2_it->second;
    
    if(v1_it->first == v2_it->first)
    {
      // (v-w)^2/(v+w) - v - w = -4 vw/(v+w)
      // we move the -4 out
      if(vi + wi != 0.0) score += vi * wi / (vi + wi);
      
      // move v1 and v2 forward
      ++v1_it;
      ++v2_it;
    }
    else if(v1_it->first < v2_it->first)
    {
      // move v1 forward
      v1_it = v1.lower_bound(v2_it->first);
    }
    else
    {
      // move v2 forward
      v2_it = v2.lower_bound(v1_it->first);
    }
  }
    
  // this takes the -4 into account
  score = 2. * score; // [0..1]

  return score;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double KLScoring::score(const BowVector &v1, const BowVector &v2) const
{ 
  BowVector::const_iterator v1_it, v2_it;
  const BowVector::const_iterator v1_end = v1.end();
  const BowVector::const_iterator v2_end = v2.end();
  
  v1_it = v1.begin();
  v2_it = v2.begin();
  
  double score = 0;
  
  // all the items or v are taken into account
  
  while(v1_it != v1_end && v2_it != v2_end)
  {
    const WordValue& vi = v1_it->second;
    const WordValue& wi = v2_it->second;
    
    if(v1_it->first == v2_it->first)
    {
      if(vi != 0 && wi != 0) score += vi * log(vi/wi);
      
      // move v1 and v2 forward
      ++v1_it;
      ++v2_it;
    }
    else if(v1_it->first < v2_it->first)
    {
      // move v1 forward
      score += vi * (log(vi) - LOG_EPS);
      ++v1_it;
    }
    else
    {
      // move v2_it forward, do not add any score
      v2_it = v2.lower_bound(v1_it->first);
      // v2_it = (first element >= v1_it.id)
    }
  }
  
  // sum rest of items of v
  for(; v1_it != v1_end; ++v1_it) 
    if(v1_it->second != 0)
      score += v1_it->second * (log(v1_it->second) - LOG_EPS);
  
  return score; // cannot be scaled
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double BhattacharyyaScoring::score(const BowVector &v1, 
  const BowVector &v2) const
{
  BowVector::const_iterator v1_it, v2_it;
  const BowVector::const_iterator v1_end = v1.end();
  const BowVector::const_iterator v2_end = v2.end();
  
  v1_it = v1.begin();
  v2_it = v2.begin();
  
  double score = 0;
  
  while(v1_it != v1_end && v2_it != v2_end)
  {
    const WordValue& vi = v1_it->second;
    const WordValue& wi = v2_it->second;
    
    if(v1_it->first == v2_it->first)
    {
      score += sqrt(vi * wi);
      
      // move v1 and v2 forward
      ++v1_it;
      ++v2_it;
    }
    else if(v1_it->first < v2_it->first)
    {
      // move v1 forward
      v1_it = v1.lower_bound(v2_it->first);
      // v1_it = (first element >= v2_it.id)
    }
    else
    {
      // move v2 forward
      v2_it = v2.lower_bound(v1_it->first);
      // v2_it = (first element >= v1_it.id)
    }
  }

  return score; // already scaled
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double DotProductScoring::score(const BowVector &v1, 
  const BowVector &v2) const
{
  BowVector::const_iterator v1_it, v2_it;
  const BowVector::const_iterator v1_end = v1.end();
  const BowVector::const_iterator v2_end = v2.end();
  
  v1_it = v1.begin();
  v2_it = v2.begin();
  
  double score = 0;
  
  while(v1_it != v1_end && v2_it != v2_end)
  {
    const WordValue& vi = v1_it->second;
    const WordValue& wi = v2_it->second;
    
    if(v1_it->first == v2_it->first)
    {
      score += vi * wi;
      
      // move v1 and v2 forward
      ++v1_it;
      ++v2_it;
    }
    else if(v1_it->first < v2_it->first)
    {
      // move v1 forward
      v1_it = v1.lower_bound(v2_it->first);
      // v1_it = (first element >= v2_it.id)
    }
    else
    {
      // move v2 forward
      v2_it = v2.lower_bound(v1_it->first);
      // v2_it = (first element >= v1_it.id)
    }
  }

  return score; // cannot scale
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

