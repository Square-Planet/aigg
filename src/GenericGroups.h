#ifndef AIGG_GENERIC_GROUPS_H
#define AIGG_GENERIC_GROUPS_H

#include <type_traits>
#include <utility>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <sstream>

#include <functional> // for std::function and std::hash
#include <algorithm> // for equal_range and stable_sort

#include "CxSequenceIo.h"


namespace groups {


// anything to ascii
template<class T>
std::string atoa(T const &a) {
   std::stringstream str;
   str << a;
   return str.str();
}


using namespace ct::sequence_io;

enum FCosetSide {
   // A left coset c X is the set c X = {c x; x ∈ X} where X ⊂ G is a
   // subgroup and c ∈ G is a coset representative
   COSET_Left,
   // A right coset X c is the set X c = {x c; x ∈ X}
   COSET_Right
};

static char const *AsString(FCosetSide Side) {
   return ((Side == COSET_Left)? "left" : ((Side == COSET_Right)? "right" : "??#unk??"));
}


// we use null_index_v<index_t> to mark invalid / yet-unassigned sequence indices of objects
template<class index_t> inline constexpr index_t null_index_v = std::is_signed_v<index_t> ? -1 : index_t(-1);

// returns whether index `i` equals null_index_v<T>
template<class index_t> bool index_assigned_q(index_t i) { return i != null_index_v<index_t>; };

// returns whether index `i` is assigned and lies in range {0, 1, ..., N-1}.
template<class index_t> bool index_valid_q(index_t i, size_t N) { return i >= 0 && size_t(i) < N && index_assigned_q(i); };
template<class index_t> bool index_invalid_q(index_t i, size_t N) { return !index_valid_q(i,N); };


// index of a group element/set element in relation an original ordered list G = (g_i)_{i=0}^{|G|-1}
using element_index_t = size_t;
// index used as identifier to the equivalence classes of cosets c X = {c x; x ∈ X}
using coset_index_t = size_t;

using element_index_list_t = std::vector<element_index_t>;
using coset_index_list_t = std::vector<coset_index_t>;


template<class FGroupElement, class FGroupElementList, class FElementKeyFn, class FGroupMultiplyFn, class FGroupInvertFn>
void _VerifyGroupStructureT(FGroupElement const &e_identity , FGroupElementList const &Elements, FGroupMultiplyFn EvalProduct, FGroupInvertFn EvalInverse, FElementKeyFn MakeElementKey);


// A helper structure used to relate elements `e` in relation
// to an initial ordered list of elements (g_i)_{i=0}^{|G|-1} of G.
//
// Basic implementation as sorted table, using binary searches
// on element keys to locate elements.
template<class FElement, class FKey>
struct FElementSearchTable {
   typedef FKey key_type;
   typedef FElement element_type;

   typedef std::function<key_type(element_type)> FMakeKeyFn;

   // Given an ordered list of elements G = (g_i; i ∈ {0,1,...}), construct a
   // search table for mapping any e ∈ G to its element index `i` (that is,
   // the integer `i` with g_i = e)
   template<class FSequence>
   explicit FElementSearchTable(FSequence const &Elements, FMakeKeyFn MakeKey);

   // If queried element `e` fulfills e = g_i for some g_i ∈ G from the
   // original list, return element index `i`. Otherwise return empty optional.
   std::optional<element_index_t> FindIndexQ(FElement const &e) const;
   // Return whether `e` is contained in original element list G.
   bool ConainsQ(FElement const &e) const { return FindIndexQ(e).has_value(); }
   // Find index `i` with e == g_i of an element `e` assumed to be contained
   // in the original sequence (g_i)_{i=0}^{|G|-1}.
   // Will raise exception if not actually there.
   element_index_t IndexQ(FElement const &e) const;

   size_t size() const { return m_SortedKiPairs.size(); }
protected:
   using FKiPair = std::pair<FKey, element_index_t>; // key, index pair
   using FSearchTable = std::vector<FKiPair>;
   static bool KiCompare_KeyOnly(FKiPair const &a, FKiPair const &b) {
      return a.first < b.first;
   }

   FSearchTable
      // maps element keys to the index of the original element in the provided
      // sequence. The array is and remains sorted.
      m_SortedKiPairs;
   FMakeKeyFn
      _KeyQ;

};


template<class FElement, class FKey>
template<class FSequence>
FElementSearchTable<FElement, FKey>::FElementSearchTable(FSequence const &Elements, FMakeKeyFn MakeKey)
   : _KeyQ(MakeKey)
{
   for (auto &&[ie,e] : enumerate(Elements))
      m_SortedKiPairs.emplace_back(_KeyQ(e), ie);
   std::stable_sort(m_SortedKiPairs.begin(), m_SortedKiPairs.end());
   for (size_t i = 1; i < m_SortedKiPairs.size(); ++ i) {
      // array is sorted now, so if there are any two equal/equivalent keys,
      // they should stand next to each other by now. So check if we found
      // any duplicates. That would be an error---they keys are supposed
      // to uniquely identify each element, and each element in the original
      // group should only come once.
      FKey const
         &ki = m_SortedKiPairs[i].first,
         &kj = m_SortedKiPairs[i-1].first;
      if (!(ki < kj || kj < ki))
         throw std::runtime_error("FElementSearchTable::c'tor(): Element list contains duplicates or sort key function is defective");
   }
}

template<class T>
struct MagicTypeRevealer;


template<class FElement, class FKey>
std::optional<element_index_t> FElementSearchTable<FElement, FKey>::FindIndexQ(FElement const &e) const
{
   FKiPair
      // compute element search key. And combine it with a null index (which is
      // just ignored by the used sorting predicate) to make it type-compatible
      // with the stored table elements
      e_and_null_index(_KeyQ(e), null_index_v<typename FKiPair::second_type>);
   auto [EqFirst, EqLast] = std::equal_range(m_SortedKiPairs.begin(), m_SortedKiPairs.end(),
      e_and_null_index, &KiCompare_KeyOnly);
   if (EqFirst == EqLast)
      // not there.
      return {};
   if (EqLast - EqFirst > 1)
      // note that the keys already were unique. If this still happens, maybe the sorting key
      // is broken (not transitive?).
      throw std::runtime_error("FElementSearchTable::FindIndexQ: input key compares equal to multiple table keys");
   FKiPair const ki = *EqFirst;
   assert(index_valid_q(ki.second, size()));
   return {ki.second};
}


template<class FElement, class FKey>
element_index_t FElementSearchTable<FElement, FKey>::IndexQ(FElement const &e) const {
   std::optional<element_index_t>
      res = FindIndexQ(e);
   if (!res)
      throw std::runtime_error("FElementSearchTable: failed to locate index of a group element which was supposed to exist.");
   element_index_t ie = res.value();
   assert(index_valid_q(ie, this->size()));
   return ie;
}




// Identify all cosets of `G`, and classify all elements g ∈ G into cosets.
//
// Returns: {iCosetReps, iElementCosets}, with:
// - iCosetReps: map of (coset index) --> (index in `gElements` of one of its representatives),
// - iElementCosets: map of (index in `gElements`) --> (assigned coset index).
//
// This function makes no attempt to assign either coset indices or coset
// representatives in any specific manner. In particular, it neither guarantees
// that the coset containing the identity element gets assigned index 0 nor that
// the identity element is chosen as one of the coset representatives.
template<class FListX, class FListG, class FTableG>
auto _AssignInitialCosets(FListX const &xElements, FCosetSide Side, FListG const &gElements,
      FTableG const &gTable, FPrintLevel PrintLevel = FPrintLevel::Off, ct::FLog *pio=0)
   -> std::tuple<element_index_list_t, coset_index_list_t>
{
   // first make the cosets themselves: we represent them with
   // a table which maps the element indices (in the original list `ParentGroupElements`)
   // to the index of the coset they belong to.
   size_t const
      Ng = gElements.size(),
      Nx = xElements.size(),
      nCosetsExpected = (assert(Ng % Nx == 0), Ng/Nx);

   coset_index_list_t
      // for each `g` element index of its assigned coset (on return),
      // or npos if not assigned yet (during incremental build up)
      iElementCosets(Ng, null_index_v<coset_index_t>);
   element_index_list_t
      // for each encountered coset, element index of one representative
      iCosetReps;
   iCosetReps.reserve(nCosetsExpected);

   // loop over representatives `c` of each coset. Each `g` element could take
   // the role, so that's where we start.
   for (element_index_t ic = 0; ic != gElements.size(); ++ ic) {
      // already assigned c to a coset?
      coset_index_t
         cCosetIndex = iElementCosets[ic];
      if (index_assigned_q(cCosetIndex)) {
         // yes, already processed in coset of earlier element
         assert(index_valid_q(cCosetIndex, nCosetsExpected));
         continue;
      }
      // `c` was not in any of the cosets generated by earlier group elements.
      // So it is a representative of a new one.
      cCosetIndex = iCosetReps.size();
      iCosetReps.push_back(ic);
      auto const &c = gElements[ic];

      // Now find all elements of the coset (i.e., all {c x; x ∈ X} or all
      // {x c; x ∈ X}), and assign them to the same coset as `c`.
      for (auto const &x : xElements) {
         typename FListG::value_type
            h = (Side == COSET_Left)? (c * x) : (x * c);
         element_index_t
            ih = gTable.IndexQ(h);
         assert(index_valid_q(ih,Ng));
         coset_index_t
            &hCosetIndex = iElementCosets[ih];
         // Confirm that hCosetIndex is still unassigned (if it were assigned,
         // then cCosetIndex should also have been assigned, as both are in the
         // same coset).
         if (index_assigned_q(hCosetIndex))
            throw std::runtime_error("_AssignInitialCosets: encountered representative"
               " of a new coset which was already assigned another coset index. X not a group?");
         hCosetIndex = cCosetIndex;
      }
      // Since `X` includes an identity element, `c` itself should not
      // need any special handling --- its index *should* have been
      // assigned in the loop just now. Make sure this is indeed the case
      cx_assert_rt(iElementCosets[ic] == cCosetIndex);
   }

   // number of cosets found should be identical to the index [G : X] of X in G
   cx_assert_rt(iCosetReps.size() == nCosetsExpected);
   // verify that all elements have indeed been assigned a coset index
   for (element_index_t ic = 0; ic != gElements.size(); ++ ic)
      cx_assert_rt(index_valid_q(iElementCosets[ic], iCosetReps.size()));

   if (PrintLevel.BasicQ() && bool(pio)) {
      pio->Write("    Sorted |G| = {} group elements {{g; g ∈ G}} into {} {} cosets of subgroup X (with |X| = {})",
         gElements.size(), iCosetReps.size(), AsString(Side), xElements.size());
   }
   return {iCosetReps, iElementCosets};
}


// Updates the element-to-coset mapping iElementCosets to be consistent with the
// ordering implied by the new coset representatives
void _UpdateCosetMapping(coset_index_list_t &iElementCosets, element_index_list_t const &iCosetReps)
{
   size_t
      Nc = iCosetReps.size();
   // Compute the mapping from the previous coset indices to the
   // new ones.
   coset_index_list_t
      iOldToNew(Nc, null_index_v<coset_index_t>);

   for (coset_index_t iNew = 0; iNew != Nc; ++ iNew) {
      // the i'th element of iCosetReps is a representative of the coset with
      // new index `i`. To found out what to update, we just need to look up
      // which coset index the new representative had originally.
      coset_index_t iOld = iElementCosets.at(iCosetReps[iNew]);
      if (index_invalid_q(iOld, Nc) || index_valid_q(iOldToNew[iOld], Nc))
         throw std::runtime_error("_UpdateCosetMapping(): encountered two new"
         " supposedly-unique coset representatives originally assigned to the same coset.");
      iOldToNew.at(iOld) = iNew;
   }
   // cross check to make sure that all previous cosets got new indices
   // assigned.
   for (coset_index_t iOld = 0; iOld != Nc; ++ iOld)
      cx_assert_rt(index_valid_q(iOldToNew[iOld], Nc));

   for (element_index_t ig = 0; ig != iElementCosets.size(); ++ ig)
      iElementCosets[ig] = iOldToNew.at(iElementCosets[ig]);
}


// Attempt to find a single generator `g` which spans a cyclic group
//
//   ⟨g⟩ = {g^i; i ∈ {0,1,...,Nc-1}}
//
// of order Nc = nCosets exactly, and with each element of ⟨g⟩ lying in a
// different coset.
// We can then use the {g^i; i ∈ {0,1,...,Nc-1}} as coset representatives,
// and the exponents `i` as new coset indices.
//
// Comments:
// - This is a very special situation, and not guaranteed to be possible at
//   all for general groups and general cosets.
// - It may be possible to generalize this function to construct reasonably
//   meaningful coset representatives also in other cases.
//
//   It may also be useful when developing an algorithm to incrementally
//   construct a set of cyclic generators for a group for which only the
//   elements are known, but not any convenient non-bruteforce way of re-
//   generating those elements. That is why I made it.
//
// Returns:
// - Empty list if no generator fulfilling the criteria was found
// - For the (g^i; i = 0,1,...,Nc-1) which span the cosets, list of
//   element indices in gElements.
template<class FListG, class FTableG, class FElementG>
element_index_list_t _TryFindCosetSpanningCyclicGenerator(coset_index_list_t &iElementCosets,
   element_index_list_t const &iCosetReps, FElementG gIdentity, FListG const &gElements, FTableG const &gTable)
{
   element_index_t
      igIdentity = gTable.IndexQ(gIdentity);

   size_t Nc = iCosetReps.size();
   std::vector<bool> bCosetCovered;
   element_index_list_t igPows; // [i]: element index in `gElements` of g^i.
//    bool Works = false;

   for (auto &&[ig,g] : enumerate(gElements)) {
      igPows.assign(Nc, null_index_v<element_index_t>);
      bCosetCovered.assign(Nc, false);
      // start search with power i=1 instead of i=0.
      // Less elegant, but avoids many unnecessary IndexQ(id) and id*g
      FElementG
         gPow = g;
      for (size_t ipow = 1; ipow < 1 + Nc; ++ ipow) {
         element_index_t igPow = gTable.IndexQ(gPow);
         if ((igPow == igIdentity) != (ipow == Nc))
            // `g` not cyclic with order `Nc` --- either
            // g^i == id already for i < Nc or g^Nc != id.
            break;
         igPows[ipow % Nc] = igPow;
         coset_index_t gPowCosetIndex = iElementCosets[igPow];
         assert(index_valid_q(gPowCosetIndex, Nc));
         if (bCosetCovered.at(gPowCosetIndex))
            // g^n is in the same coset as an earlier g^i. This is not
            // the generator we are looking for.
            break;
         bCosetCovered[gPowCosetIndex] = true;
         if (ipow == Nc) {
            assert(std::all_of(bCosetCovered.begin(), bCosetCovered.end(), [](bool b) {return b;}));
            assert(std::all_of(igPows.begin(), igPows.end(), [Nc](bool ig) {return ig < Nc;}));
            // if we got to this point, then g is a cyclic generator
            // of order Nc and all its powers lie in different cosets.
            return igPows;
         }
         gPow = gPow * g;
      }
   }
   return {};
}


static bool _IsPrimeNumber(size_t n) {
   if (n < 2) return false;
   // amazingly clever test of primeness, right?
   for (size_t i = 2; i < n; ++ i) {
      if (n % i == 0)
         return false;
   }
   return true;
}


// Find the single generator `g` which spans largest order cyclic group
//
//   ⟨g⟩ = {g^i; i ∈ {0,1,...,Nc-1}}
//
// with each element of ⟨g⟩ lying in a different coset of X in G, and
//
//   (~g * x * g) ∈ X
//
// In that case, ⟨g + X⟩ is also a group, and each of its elements
// can be written as a product of x ∈ X and one of the coset representatives
// {g^i; i ∈ {0,1,...,Nc-1}}
//
// Returns:
// - Empty list if no generator fulfilling the criteria was found
// - For the (g^i; i = 0,1,...,Nc-1) which span the cosets, list of
//   element indices in gElements.
template<class FListG, class FTableG, class FListX, class FElementG>
element_index_list_t _FindLargestCyclicGeneratorGroupExtension(coset_index_list_t &iElementCosets,
   element_index_list_t const &iCosetReps, FElementG gIdentity, FListG const &gElements,
   FTableG const &gTable, FListX const &xElements,
   FPrintLevel PrintLevel = FPrintLevel::Off, ct::FLog *pio=0)
{
   element_index_t
      igIdentity = gTable.IndexQ(gIdentity);
   coset_index_t
      idCosetIndex = iElementCosets.at(igIdentity);

   size_t Nc = iCosetReps.size();
   std::vector<bool> bCosetCovered;
   element_index_list_t igPows; // [i]: element index in `gElements` of g^i.
//    bool Works = false;

   // element indices of {g^i; i = 0,1,...} for the cyclic generator
   // covering the largest number of distinct cosets so far.
   element_index_list_t igPows_BestSoFar;

   for (auto &&[ig,g] : enumerate(gElements)) {
      igPows.clear();
      igPows.push_back(igIdentity);
      bCosetCovered.assign(Nc, false);
      bCosetCovered[idCosetIndex] = true;
      // Compute g^i for i = 1, 2, ... (skipping i = 0, as we already added id).
      // For as long as all lie in different cosets (--> at most Nc loop iterations)
      FElementG
         gPow = g;
      for (size_t ipow = 1; ; ++ ipow) {
         element_index_t
            igPow = gTable.IndexQ(gPow);
         coset_index_t
            gPowCosetIndex = iElementCosets[igPow];
         assert(index_valid_q(gPowCosetIndex, Nc));
         assert(index_valid_q(gPowCosetIndex, iElementCosets.size()));

         if (bCosetCovered.at(gPowCosetIndex)) {
            // g^n is in the same coset as an earlier g^i.

            // Since we are looking for cyclic generators, we only want to
            // consider this one if the first time it happens is when we hit
            // back the identity element again.
            if (igPow != igIdentity)
               // ...not what we're looking for.
               igPows.clear();
            break;
         }
         igPows.push_back(igPow);
         bCosetCovered[gPowCosetIndex] = true;
         gPow = gPow * g;
      }
      // can't have generator images in more distinct cosets than there are cosets.
      assert(igPows.size() <= Nc);
      if (igPows.size() == Nc) {
         // if the generator orbit has Nc entries in distinct cosets, then it
         // spans *all* the cosets of X in G, and therefore, combined with X,
         // the entire group. This implies the conjugation property checked
         // below without testing.
         igPows_BestSoFar.swap(igPows);
         break;
      }

      if (false && !_IsPrimeNumber(igPows.size()))
         continue;

      if (igPows.size() > igPows_BestSoFar.size()) {
         // still needs checking: if all ~g*x*g ∈ X
         // That is the case if ~g*x*g lies in the coset containing the identity
         // element (which is identical to X itself)
         for (auto const &x : xElements) {
            // ^-- one also could get the x elements by iterating over g ∈ G and
            // skipping over g with coset indices differing from id's. However,
            // we may use different data types for X and G (e.g., permuflections
            // and rotation matrices), and it might be useful to retain this ability.
            FElementG
               xConj = ~g * x * g;
            coset_index_t
               xConjCosetIndex = iElementCosets[gTable.IndexQ(xConj)];
            if (xConjCosetIndex != idCosetIndex) {
               igPows.clear();
               break;
            }
         }
      }
      if (igPows.size() > igPows_BestSoFar.size()) {
         // At this moment, igPows should contain:
         // - element indices of {g^i; i = 0,1,...,m-1},
         // - all lying in different cosets,
         // - and g^m would be identity if included (and therefore
         //   not be in a distinct coset anymore, since we have g^0 == id)
         // If it has more elements than the last such set, remember it.
         if (igPows.size() > igPows_BestSoFar.size())
            igPows_BestSoFar.swap(igPows);
      }

   }
   if (PrintLevel.BasicQ()) {
      if (igPows_BestSoFar.empty()) {
         pio->Write("    Failed to find cyclic generator extending |X| = {} to a larger subgroup of |G| = {}",
            xElements.size(), gElements.size());
      } else {
         pio->Write("    Found cyclic generator `c ∈ (G ∖ X)` with |⟨c⟩| = {}, c X = X c, and all c^i (i=0,...,{}) in distinct cosets (of [G:X] = {})",
            igPows_BestSoFar.size(), igPows_BestSoFar.size()-1, iCosetReps.size());
      }
   }

   // note: this may still be empty. For a group consisting of
   // identity only, it should consist of only the id element.
   return igPows_BestSoFar;
}



template<class FListG, class FTableG, class FElementG>
void _VerifyGroupStructure1(FElementG gIdentity, FListG const &gElements, FTableG const &gTable)
{
   using G = FElementG;
   _VerifyGroupStructureT(gIdentity, gElements,
      [](G const &a, G const b){ return a * b; },
      [](G const &a){ return ~a; },
      [&gTable](G const &a){ return gTable.IndexQ(a); }
   );
}


// index of generator group element in parent list, and cyclic order
// (the smallest n > 0 with g^n == id)
typedef std::pair<element_index_t, size_t>
   generator_info_t;

template<class FListG, class FTableG, class FElementG>
auto _TryBuildCyclicGeneratorChain(FElementG gIdentity, FListG const &gElements,
      FTableG const &gTable, FPrintLevel PrintLevel = FPrintLevel::Off, ct::FLog *pio=0)
   -> std::vector<generator_info_t>
{
   if (PrintLevel.MoreQ())
      pio->Write("\n ░░ NEXT: Search cyclic generator chain for group G with {} elements", gElements.size());
//    element_index_t
//       igIdentity = gTable.IndexQ(gIdentity);
   // start with the identity group. Yes, not very imaginative.
   FListG
      xElements = {gIdentity};
   std::vector<generator_info_t>
      CyclicGens;

   while (xElements.size() < gElements.size()) {
      // construct cosets for current X in G
      auto [iCosetReps, iElementCosets] = groups::_AssignInitialCosets(
         xElements, COSET_Left, gElements, gTable, PrintLevel-1, pio);

      // try find a cyclic generator h of order > 1 for which {h x ∈ X} == ⟨h + x ∈ X⟩
      element_index_list_t
         ihPows = _FindLargestCyclicGeneratorGroupExtension(
            iElementCosets, iCosetReps, gIdentity, gElements, gTable, xElements, PrintLevel-1, pio);
      // ^-- returns element indices in gTable for h^i
      if (ihPows.size() < 2) {
         // construction failed --- there are still some group elements
         // left, and we did not manage to find a cyclic generator
         // with X ⊂ h X == ⟨X + h⟩ ⊆ G
         pio->Write("    FAILED: Unable to extend subgroup X ⊂ G with |X| = {} to full G with |G| = {} elements. Could not find additional cyclic generator", xElements.size(), gElements.size());
         return {};
      }
      // append h^1 = h to generator list; the function returns successfully, if
      // for n == ihPows.size(), we have h^n == e and all lower h^i are in
      // distinct cosets.
      CyclicGens.push_back({ihPows[1], ihPows.size()});
      // If a `h` fulfilling all previous conditions was found, then with n is
      // the generator order. We can make the next larger group simply as
      //
      //     Y := {h^i; x ∈ X, i ∈ {0,1,...,(n-1)}}
      FListG
         yElements;
      for (auto const &x : xElements)
         for (element_index_t ihpow : ihPows)
            yElements.push_back(gElements[ihpow] * x);

      if (PrintLevel.MoreQ() || (PrintLevel.BasicQ() && yElements.size() == gElements.size())) {
         std::stringstream strg, stre;
         for (auto &&[igen,ih_hord] : enumerate(CyclicGens)) {
            auto &&[ih,hord] = ih_hord; // hm... no multi-level structured binding :(
            if (igen != 0) {
               strg << " "; stre << ", ";
            }
            char idx = char('i' + igen);
            strg << "h^" << idx;
            stre << idx << "∈{0,…," << (hord-1) << "}";
//             stre << idx << "∈{1,…," << hord << "}";
            // ^-- I'd probably better span it as 0...(n-1), but 1...n should also work,
            //     since h^n == id, too
         }

         pio->Write("    Constructed subgroup X ⊂ G with |X| = {} of |G| = {} elements, generated as X = {{{}; {}}} [FIXME: reversed!]",
            yElements.size(), gElements.size(), strg.str(), stre.str());
         if (PrintLevel.MoreQ())
            pio->WriteLine();
      }

#ifdef _DEBUG
      _VerifyGroupStructure1(gIdentity, yElements, gTable);
#endif // _DEBUG

      // use Y as the next subgroup X
      xElements.swap(yElements);
   }
   assert(xElements.size() == gElements.size());
   return CyclicGens;
}









namespace _VerifyGroupStructureHelper {
   static void Fail(std::string const &reason) {
      throw std::logic_error("_VerifyGroupStructureT(): FAILED: " + reason);
   }

   template<class FGroupElement, class... Args>
   static void Fail(std::string const &reason, FGroupElement const &e, Args ...args) {
      return Fail(reason + "\n [-"+atoa(1+sizeof...(Args))+"] Offending element e: " + atoa(e), (args)...);
   }
   // ^--- apparently not allowed in either local class or local function :/

   template <class T> T &&return_argument_unchanged(T &&obj) { return std::forward(obj); }
}; // std::identity<FGroupElement> // <-- >= C++20 only.
//  + " //\t" + atoa(MakeSortKey_RotationMatrix<double>(e))


// Explicitly tests whether all axioms which a group are supposed to fulfill
// are indeed fulfilled for a given (finite) set of group elements and group multiplication.
// Also tests if the supplied identity element, inverse functions, and element -> key
// mapping functions behave as they should.
//
// Arguments:
// - e_identity: identity element of the group (i.e., element with ∀g: eid * g == g)
// - Elements: full set of {g; g ∈ G} which are believed to be all the group elements
// - EvalProduct: function object callable with two arguments g,h ∈ G, and returning
//   the product (g*h).
// - EvalInverse: function object callable with one argument g ∈ G and returning the
//   inverse element g⁻¹
// - MakeElementKey: function object callable with one argument g ∈ G, and
//   returning either `g` itself, or another object to be used in g's stead
//   as a well-behaved set key (in particular, the key objects must be sortable and
//   define a strict strong ordering on the group elements;
//   For all group elements g,h ∈ G it must fulfill: (key(g) == key(h)) == (g == h)
template<class FGroupElement, class FGroupElementList, class FElementKeyFn, class FGroupMultiplyFn, class FGroupInvertFn>
void _VerifyGroupStructureT(FGroupElement const &e_identity , FGroupElementList const &Elements, FGroupMultiplyFn EvalProduct, FGroupInvertFn EvalInverse, FElementKeyFn MakeElementKey)
{
   namespace loc = _VerifyGroupStructureHelper;
   if (Elements.empty())
      return;
   static_assert(std::is_invocable_v<decltype(MakeElementKey), decltype(*Elements.begin())>,
      "_VerifyGroupStructureT: Required set key generation function of group elements cannot be invoked for provided group elements");

   // collect keys of group elements into a set (to simplify lookup of existence of elements)
   typedef std::remove_reference_t<decltype(MakeElementKey(*Elements.begin()))>
      FElementKey;
   std::set<FElementKey>
      ElementKeys;
   for (FGroupElement const &e : Elements)
      ElementKeys.insert(MakeElementKey(e));
   if (ElementKeys.size() != Elements.size())
      loc::Fail("Element list contains duplicates or sort key function is defective");

   auto ContainsQ = [&](FGroupElement const &e) {
      auto const &e_key = MakeElementKey(e);
      return ElementKeys.count(e_key) == 1;
   };
   auto EqualQ = [&](FGroupElement const &e, FGroupElement const &f) {
      auto const &e_key = MakeElementKey(e);
      auto const &f_key = MakeElementKey(f);
      return e_key == f_key;
//       return !(e_key < f_key || f_key < e_key);
   };


   // checks regarding identity element
   if (!ContainsQ(e_identity))
      loc::Fail("identity element not in group", e_identity);
   for (FGroupElement const &e : Elements) {
      // check if identity element behaves as ∀e ∈ G:  e * id == e
      FGroupElement also_e = EvalProduct(e, e_identity);
      if (!EqualQ(e, also_e))
         loc::Fail("provided identity element does not fulfill e * id == e", e);
   }

   // checks regarding inverse elements
   for (FGroupElement const &e : Elements) {
      FGroupElement inv_e = EvalInverse(e);
      if (!ContainsQ(inv_e))
         loc::Fail("inverse element e^{-1} not in group", e);
      if (!EqualQ(EvalProduct(e, inv_e), e_identity))
         loc::Fail("inverse element defective: e * e^{-1} != id", e);
   }

   // check if group is closed under group multiplication
   for (FGroupElement const &e : Elements)
      for (FGroupElement const &f : Elements)
         if (!ContainsQ(EvalProduct(e,f)))
            loc::Fail("group not closed: e * f not in group", e, f);

   // check if group multiplication is associative (e f) h == e (f h)
   for (FGroupElement const &e : Elements)
      for (FGroupElement const &f : Elements)
         for (FGroupElement const &g : Elements) {
            FGroupElement
               e_fg = EvalProduct(e, EvalProduct(f, g)),
               ef_g = EvalProduct(EvalProduct(e, f), g);
            if (!EqualQ(e_fg, ef_g))
               loc::Fail("group multiplication not associative: (e*f)*g != e*(f*g)", e, f, g); // , e_fg, ef_g
         }
//    io.Write("_VerifyGroupStructureT(): {} element group passed. [`{}`]", Elements.size(), typeid(e_identity).name());
//    throw std::runtime_error("absicht.");
}



} // namespace groups


#endif // AIGG_GENERIC_GROUPS_H
