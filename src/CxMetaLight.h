/// @file CxMetaLight.h C++ template/parameter pack meta-programming primitives
/// (light-weight parts only, i.e., no #includes of `<tuple>`, `<array>`,
/// `<functional>`, etc); see `CxMeta.h` for overview and documentation.
#ifndef CX_META_LIGHT_H
#define CX_META_LIGHT_H

#include <type_traits>
#include <iterator> // for defining iterator categories: https://en.cppreference.com/w/cpp/iterator/iterator_tags
#include <utility> // for std::move/std::forward; this should also define std::pair, and therefore std::tuple_size, std::tuple_element, and std::get.


// #include <functional> // for reference wrapper (more precisely, the associated unwrapping with reduced decay)


// namespace for the c++ meta programming interface
#define _CxMeta cx
#define _Namespace_CxMeta_Begin namespace cx
#define _Namespace_CxMeta_End

// ...and add some template-programming stuff, too. One can argue that this is
// indeed about iteration, although the nature is a bit different, of course.
//
// Note: This is intentionally put into a separate namespace (not ct). I have
// often used completely different things called "index lists"....
_Namespace_CxMeta_Begin {
   // Make an implementation of something like std::index_sequence, which, unfortunately,
   // is >= C++14 only. E.g., MakeIndexList<5>() yields a return type
   // ct::TIndexList<0, 1, 2, 3, 4>, which can be used in template<class
   // TIndexList<Indices...> > in variable template argument lists or similar.
   //
   // Note to self: A common usage pattern of this works like that:
   //
   //       template<unsigned... Is, class FTupleLike>
   //       void WriteTupleImpl(std::ostream &out, FTupleLike const &A, index_list_t<Is...>) const {
   //          using std::get;
   //          int dummy[] = {
   //             ((out << ((Is == 0)? "" : get<1>(m_SequenceSeptors)) << get<Is>(A)), 0)...
   //             // ^-- the (stuff,0)... which is repeated is the array element, each of
   //             //     which evaluates to 0 because this is a comma-operator. The array
   //             //     thing is apparently the canonical way in C++11 to execute a
   //             //     statment in input order for each entry in the parameter pack.
   //          };
   //          (void)dummy; // suppress unused warning.
   //       }
   //
   //       template<class... Args>
   //       std::ostream &WriteTuple(std::ostream &out, std::tuple<Args...> const &A) const {
   //          // ...
   //          WriteTupleImpl(out, A, index_list_for_t<Args...>());
   //          // ...
   //       }
   //
   // There are two functions: The originally called one (here: WriteTuple) just
   // creates the template index list, and then defers to a second function
   // which takes it as parameter pack template argument (here: WriteTupleImpl).
   // Inside this function, one can then pack-expand the index list template
   // (here: Is) to get a comma- separate sequence for all indices.
   //
   // E.g. for Is = <0, 1, 2> (index list of three compile-time arguments), we'd
   // see those expansions:
   //
   //    fragment `m[Is]...`  --> turns into `m[0], m[1], m[2]`
   //    fragment `func(m[Is]...)` -->  turns into  `func(m[0], m[1], m[2])`
   //    fragment `func((m[Is])...)` -->  turns into  `func((m[0]), (m[1]), (m[2]))`
   //    // ...i.e., (doesn't do anything here, but sometimes the extra parenthesis are needed to disambiguate)
   //    fragment `(out << get<Is>(a))...` --> turns into `(out << get<0>(a)), (out << get<1>(a)), (out << get<2>(a))`
   //    fragment `(func(x,Is,y), 0)...` --> turns into `(func(x,0,y), 0), (func(x,1,y), 0), (func(x,2,y), 0)`
   //
   // etc. So the pack expansions are similar to text replacements: take
   // everything in the expression right before the "..." (use parenthesis to
   // define enclosed scope in detail) and call it "Expr". Then for each entry
   // in the index list expansion, one gets a new entry in which the occurrence
   // of the index list (here `Is`) in the expression is replaced by its value,
   // and these new expressions are written next to each other with a comma
   // between any two---regardless of context.
   //

   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: Generation of index lists as template arguments:
   // ▒▒         index_list_t<index_t...>,
   // ▒▒         index_list_for_t<class... Args>(),
   // ▒▒         index_list_of_size_t<size_t N>,
   // ▒▒         index_list_for_range_t<size_t iFirst, size_t iLast>
   // ──────────────────────────────────────────────────────────────────────────

   /// Element of template-index lists (not itself templated to reduce bloat)
   typedef unsigned index_t;
   namespace iter_detail {
      template<index_t... Integers>
      struct TIndexList {};

      template<index_t iFirst, index_t N, index_t... OtherArgs>
      struct _TIndexListHelper {
         static_assert(N <= 200000, "probably an overflow...");
         static_assert(N != 0, "should not be here if N = 0");
         typedef typename _TIndexListHelper<iFirst, (N-1), (iFirst+N-1), OtherArgs...>::type type;
      };

      template<index_t iFirst, index_t... OtherArgs>
      struct _TIndexListHelper<iFirst, 0, OtherArgs...> {
         typedef TIndexList<OtherArgs...> type;
      };
   }


   /// Primary index list type for interfaces. It is a proxy for `cx::iter_detail::TIndexList<Args...>`
   template<index_t... Args>
   using index_list_t = iter_detail::TIndexList<Args...>;

   /// Returns number of elements in index list
   template<index_t... Args>
   constexpr index_t _IndexCountQ(index_list_t<Args...>) { return sizeof...(Args); }

   /// Type evaluates to `cx::index_list_of_size_t<sizeof...(Args)>`
   template<class... Args>
   using index_list_for_t = typename iter_detail::_TIndexListHelper<0, sizeof...(Args)>::type;
   // ^-- this typedef is sometimes easier to use than the MakeIndexList
   //     alternative. A std::index_sequence_for is also in std C++, but only
   //     starting at C++17.
   template<size_t N>
   using index_list_of_size_t = typename iter_detail::_TIndexListHelper<0,N>::type;

   /// If `iFirst < iLast`, type evaluates to `index_list_t<iFirst, iFirst+1, ..., iLast-1>`,
   /// otherwise an empty integer sequence `index_list_t<>`.
   template<size_t iFirst, size_t iLast>
   using index_list_for_range_t = typename iter_detail::_TIndexListHelper<iFirst, ((iLast > iFirst) ? (iLast - iFirst) : 0)>::type;

#ifdef INCLUDE_ABANDONED
//    template<index_t N>
//    typename iter_detail::_TIndexListHelper<0,N>::type make_index_list() {
// //       return typename iter_detail::_TIndexListHelper<0,N>::type();
//       return index_list_of_size_t<N>();
//    }
//    // ^-- that still used anywhere?
#endif // INCLUDE_ABANDONED
} _Namespace_CxMeta_End


_Namespace_CxMeta_Begin {
   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: more index list stuff. List element test, and functions for various
   // ▒▒       decompositions / modifications of the lists. Not super elegant, but
   // ▒▒       all this stuff also works in C++11.
   // ──────────────────────────────────────────────────────────────────────────

   typedef unsigned index_t;

   namespace iter_detail {
      /// Provides a constexpr function
      /// ``
      ///    template _IsEqualToAny<i0, i1, i2, ...>::eval(j)
      /// ``
      /// which returns whether the function argument `j` compares equal (e.g., j == i0)
      /// to any of the template arguments i0, i1, ...
      template<index_t...> struct _IsEqualToAny;
      template<> struct _IsEqualToAny<> {
         constexpr static bool eval(index_t) { return false; }
      };
      template<index_t ref, index_t... other_ref> struct _IsEqualToAny<ref, other_ref...> {
         constexpr static bool eval(index_t trial) { return (trial == ref) || _IsEqualToAny<other_ref...>::eval(trial); }
      };
   }

   /// `ContainsQ<Is...>(trial)` returns whether any of the integers in parameter pack `Is...` is equal to `trial`
   template<index_t... Integers>
   constexpr bool ContainsQ(index_t trial) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }

   /// `ContainsQ(trial, index_list)` returns whether any of the integers in `cx::index_list_t<...> index_list` is equal to `trial`
   template<index_t... Integers>
   constexpr bool ContainsQ(index_t trial, index_list_t<Integers...>) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }

   /// `ContainsQ<IndexList>(trial)` returns whether any of the integers in `IndexList = cx::index_list_t<...>` is equal to `trial`
   template<class FIndexList>
   constexpr bool ContainsQ(index_t trial) { return ContainsQ(trial, FIndexList()); }

   namespace iter_detail {
      /// Provides a constexpr function
      /// ``
      ///    template _GetArgN<iArg>::eval(a0, a1, a2, ...)
      /// ``
      /// which returns the iArg'th function argument. Will yield compiler error
      /// if either there are no function arguments, or less function arguments than needed for iArg.
      ///
      /// Note: no r-value refs, forwarding etc. Arguments are returned by value
      /// (meant to be used with integer lists etc), to allow constexpr
      /// evaluation (even in C++11)
      template<index_t iEntry> struct _GetArgN {
         using _GetArgN_Next = _GetArgN<iEntry-1>;
         template<class Arg0, class... OtherArgs>
         constexpr static auto eval(Arg0, OtherArgs... other_args) -> decltype(_GetArgN_Next::eval(other_args...)) {
            return _GetArgN_Next::eval(other_args...);
         }
      };
      template<> struct _GetArgN<0> {
         template<class Arg0, class... OtherArgs>
         constexpr static Arg0 eval(Arg0 arg0, OtherArgs...) { return arg0; }
      };
   }

   /// Compile-time index list with basic querying functionality (C++11-compatible) to
   /// facilitate derivative constructions
   template<index_t... Integers>
   struct TIndexListQf {
      /// Number of controlled integer template arguments
      enum { nEntries = sizeof...(Integers) };
      /// Returns whether `trial` is an element of this list (i.e., coincides
      /// with any of the indices in the template argument list)
      constexpr static bool ContainsQ(index_t trial) { return iter_detail::_IsEqualToAny<Integers...>::eval(trial); }
      /// Returns `iArg`'th element of the list
      template<index_t iArg> constexpr static index_t GetArgN() { return iter_detail::_GetArgN<iArg>::eval(Integers...); }
   };

   /// convert a basic TIndexList to a TIndexListQf with query functionality (variable version)
   template<index_t... Integers>
   TIndexListQf<Integers...> _ConvertQfil(iter_detail::TIndexList<Integers...>) { return TIndexListQf<Integers...>{}; }
   /// convert a TIndexListQf with query functionality to a basic TIndexList (with no functionality) (variable version)
   template<index_t... Integers>
   iter_detail::TIndexList<Integers...> _ConvertNfil(TIndexListQf<Integers...>) { return iter_detail::TIndexList<Integers...>{}; }
   /// convert a basic TIndexList to a TIndexListQf with query functionality (type version)
   template<class FIndexList1> using _convert_qfil_t = decltype(_ConvertQfil(std::declval<FIndexList1>()));
   /// convert a TIndexListQf with query functionality to a basic TIndexList (type version)
   template<class FIndexList1> using _convert_nfil_t = decltype(_ConvertNfil(std::declval<FIndexList1>()));

   namespace iter_detail {
      using std::get;
      template<class BaseIndices, class RetainQ, index_t N = BaseIndices::nEntries, index_t... OtherArgs>
//       template<class BaseIndices, class RetainQ, index_t N = _IndexCountQ(BaseIndices()), index_t... OtherArgs>
      struct _FilterIndices {
//          static constexpr index_t iBaseIndex = BaseIndices::template GetArgN<N-1>();
         enum { iBaseIndex = BaseIndices::template GetArgN<N-1>() };
         using type = typename std::conditional<
            RetainQ::eval(N-1, iBaseIndex),
            typename _FilterIndices <BaseIndices, RetainQ, (N-1), index_t(iBaseIndex), OtherArgs...>::type,
            typename _FilterIndices <BaseIndices, RetainQ, (N-1), OtherArgs...>::type
         >::type;
      };

      template<class BaseIndices, class RetainQ, index_t... OtherArgs>
      struct _FilterIndices <BaseIndices, RetainQ, 0, OtherArgs...> {
//          typedef TIndexListQf<OtherArgs...> type;
         typedef iter_detail::TIndexList<OtherArgs...> type;
      };

      template<class Indices> struct _IsContainedIn { constexpr static bool eval(index_t, index_t value) { return Indices::ContainsQ(value); } };
      template<class Indices> struct _IsDisjunctFrom { constexpr static bool eval(index_t, index_t value) { return !Indices::ContainsQ(value); } };
      template<index_t SlotN> struct _TakeFirstN { constexpr static bool eval(index_t pos, index_t) { return pos < SlotN; } };
      template<index_t SlotN> struct _DropFirstN { constexpr static bool eval(index_t pos, index_t) { return pos >= SlotN; } };
      // ^-- could also do things like testing for index values instead of slot indices. not sure if needed atm.
   }

   /// `index_list_reject2nd_t<B = TIndexListQf<...>, E = TIndexListQf<...> >`
   /// evaluates to another TIndexListQf index list which contains all indices
   /// of `B`, in original order, provided they are *NOT* included in `E`
   template<class BaseIndices, class Indices2>
   using index_list_reject2nd_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_IsDisjunctFrom<_convert_qfil_t<Indices2> > >::type;

   /// `index_list_retain2nd_t<B = TIndexListQf<...>, I = TIndexListQf<...> >`
   /// evaluates to another TIndexListQf index list which contains all indices
   /// of `B`, in original order, provided they are *ALSO* included in `I`
   template<class BaseIndices, class Indices2>
   using index_list_retain2nd_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_IsContainedIn<_convert_qfil_t<Indices2> > >::type;

   /// Retain the first `N` indices in `BaseIndices`.
   /// @param BaseIndices A TIndexListQf<...> object with template arguments
   ///        providing the actual index list
   /// @param N maximum number of indices to retain
   /// @returns TIndexListQf of length <= N containing up to N of the first
   ///        indices in BaseIndices
   template<class BaseIndices, index_t N>
   using index_list_take_first_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_TakeFirstN<N> >::type;

   /// Evaluates the index list remaining after dropping the first `N` indices of `BaseIndices`.
   /// `BaseIndices` and return value are TIndexListQf templates.
   template<class BaseIndices, index_t N>
   using index_list_drop_first_t = typename iter_detail::_FilterIndices<_convert_qfil_t<BaseIndices>, iter_detail::_DropFirstN<N> >::type;

} _Namespace_CxMeta_End


_Namespace_CxMeta_Begin {
#ifdef INCLUDE_ABANDONED
   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: index-based extraction of values and types in parameter packs
   // ▒▒       (nth_type_of, nth_value_of)
   // ──────────────────────────────────────────────────────────────────────────

//    /// extract n'th element type of a parameter pack
//    template<int N, class... Ts>
//    using nth_type_of =
//       typename std::tuple_element<N, std::tuple<Ts...> >::type;
//
//    /// nth_value_of<I>(a0, a1,...) returns the I'th of its function arguments.
//    /// (useful for extracting n'th element value of a parameter pack)
//    /// Note: might be better of as decltype(auto) in >= C++14.
//    template<int N, class... Ts>
//    auto nth_value_of(Ts&&... ts)
//    -> nth_type_of<N, Ts...> {
//       using std::get; return get<N>(std::forward_as_tuple(ts...));
//    }
   // ^-- moved to CxMeta.h
#endif // INCLUDE_ABANDONED


   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: reductions over paramer packs (fold_left_pp / fold_right_pp,
   // ▒▒       fold_left_ipp, fold_right_ipp, sum_pp, any_pp, all_pp, etc.)
   // ──────────────────────────────────────────────────────────────────────────

   template<class result_t, class FBinaryOp>
   result_t constexpr fold_left_pp(FBinaryOp /*op*/, result_t accu) { return accu; }

   /// fold_left(op, neutral, [a,b,c,...]) := (...(((((neutral ⋇ a) ⋇ b) ⋇ c) ⋇ d) ⋇ e) ⋇ ...)
   ///
   /// Notes: _pp -> parameter pack. Output type result_t is derived from first argument.
   /// (this first argument is also what you will get if there are no `real` arguments)
   template<class result_t, class FBinaryOp, class Arg0, class... Args>
   result_t constexpr fold_left_pp(FBinaryOp op, result_t accu, Arg0 &&a0, Args&&... args) {
      return fold_left_pp<result_t, FBinaryOp>(
            op,
            op(std::move(accu), std::forward<Arg0>(a0)),
            std::forward<Args>(args)...
      );
   }


   template<class result_t, class FBinaryOp>
   result_t constexpr fold_right_pp(FBinaryOp /*op*/, result_t neutral) { return neutral; }

   /// fold_right(op, neutral, [a,b,c,...]) := (a ⋇ (b ⋇ (c ⋇ (d ⋇ (e ⋇ (...(... ⋇ neutral)...))))))
   template<class result_t, class FBinaryOp, class Arg0, class... Args>
   result_t constexpr fold_right_pp(FBinaryOp op, result_t neutral, Arg0 &&a0, Args&&... args) {
      return op(
         std::forward<Arg0>(a0),
         fold_right_pp<result_t, FBinaryOp>(op, std::move(neutral), std::forward<Args>(args)...)
      );
   }

   namespace iter_detail {
      template<class result_t, class FBinaryOp, unsigned iArg0>
      result_t constexpr _FoldLeftWithIndexHelper(FBinaryOp /*op*/, result_t accu) {
         return accu;
      }

      template<class result_t, class FBinaryOp, unsigned iArg0, class Arg0, class... Args>
      result_t constexpr _FoldLeftWithIndexHelper(FBinaryOp op, result_t accu, Arg0 &&a0, Args&&... args) {
         return _FoldLeftWithIndexHelper<result_t, FBinaryOp, iArg0+1>(
               op,
               op.template apply<iArg0>(std::move(accu), std::forward<Arg0>(a0)),
               std::forward<Args...>(args)...
         );
      }
   }  // ^-- both the std::move and the missing forwards look weird. should be re-checked.
   // _ipp: indexed parameter pack
   template<class result_t, class FBinaryOp, class... Args>
   result_t constexpr fold_left_ipp(FBinaryOp op, result_t neutral, Args&&... args) {
      return iter_detail::_FoldLeftWithIndexHelper<result_t,FBinaryOp,0u>(op, std::move(neutral), std::forward<Args>(args)...);
   }

#ifdef INCLUDE_ABANDONED
//    namespace iter_detail {
//       struct plus_ce {
//          template<class Ta, class Tb>
//          auto constexpr operator()(Ta &&a, Tb &&b) const -> decltype(a + b) { return std::forward<Ta>(a) + std::forward<Tb>(b); }
//       };
//    }
//
//    template<class T, class... Args>
//    T constexpr sum_pp(T neutral, Args&&... args) {
// //       return fold_left_pp<T>(std::plus<T>(), neutral, args...);
//       // ^-- hmpf. std::plus's operator() is not constexpr until C++14.
//       return fold_left_pp<T>(iter_detail::plus_ce(), std::move(neutral), std::forward<Args>(args)...);
//    }
//    // ^-- this works fine, and is more flexible than the simple one below in that it
//    //     can propagate return types, rather than just assuming one. However,
//    //     for most applications it might be a bit of overkill.
#endif // INCLUDE_ABANDONED

   bool constexpr any_pp() { return false; }
   bool constexpr all_pp() { return true; }
   /// any_pp(v0, v1, ...) returns if `bool(vi) == true` holds for at least one of its arguments
   template<class... Args> bool constexpr any_pp(bool a0, Args&&... args) { return a0 || any_pp(bool(args)...); }
   /// all_pp(v0, v1, ...) returns if `bool(vi) == true` holds for each of its arguments
   template<class... Args> bool constexpr all_pp(bool a0, Args&&... args) { return a0 && all_pp(bool(args)...); }

   template<class T> T constexpr sum_pp(T neutral) { return neutral; }
   template<class T, class Arg0, class... Args> T constexpr sum_pp(T neutral, Arg0 &&arg0, Args&&...args) {
      return sum_pp(std::move(neutral) + std::forward<Arg0>(arg0), std::forward<Args>(args)...);
   }

   template<class T> T constexpr prod_pp(T neutral) { return neutral; }
   template<class T, class Arg0, class... Args> T constexpr prod_pp(T neutral, Arg0 &&arg0, Args&&...args) {
      return prod_pp(std::move(neutral) + std::forward<Arg0>(arg0), std::forward<Args>(args)...);
   }

#ifdef INCLUDE_ABANDONED
//    // csum_pp = coefficient sum. point of this is that (a) all summand ars are
//    // the same type, and (b) that type is a simple scalar type which is best
//    // passed directly, by value.
//    template<class T> T constexpr csum_pp() { return T(0); }
//    template<class T> T constexpr csum_pp(T a0) { return a0; }
//    template<class T, class... Ts> T constexpr csum_pp(T a0, T a1, Ts... as) { return csum_pp(a0 + a1, as...); }
//
//    template<class T> T constexpr cprod_pp() { return T(1); }
//    template<class T> T constexpr cprod_pp(T a0) { return a0; }
//    template<class T, class... Ts> T constexpr cprod_pp(T a0, T a1, Ts... as) { return cprod_pp(a0 * a1, as...); }
#endif // INCLUDE_ABANDONED
   // csum_pp = coefficient sum. point of this is that (a) all summand ars are
   // the same type, and (b) that type is a simple scalar type which is best
   // passed directly, by value.
   template<class T> T constexpr csum_pp() { return T(0); }
   template<class T, class... Ts> T constexpr csum_pp(T a0, Ts... as) { return a0 + csum_pp<T>(as...); }

   template<class T> T constexpr cprod_pp() { return T(1); }
   template<class T, class... Ts> T constexpr cprod_pp(T a0, Ts... as) { return a0 * cprod_pp<T>(as...); }


   template<class T> T constexpr _MinQ(T const &a, T const &b) { return (a < b)? a : b; }
   template<class T> T constexpr _MaxQ(T const &a, T const &b) { return (b < a)? a : b; }
   template<class T> T constexpr min_pp(T const &a0) { return a0; }
   template<class T, class... Args> T constexpr min_pp(T const &a0, T const &a1, Args const&... args) { return _MinQ<T>(a0, min_pp(a1, args...)); }
   template<class T> T constexpr max_pp(T const &a0) { return a0; }
   template<class T, class... Args> T constexpr max_pp(T const &a0, T const &a1, Args const&... args) { return _MaxQ<T>(a0, max_pp(a1, args...)); }

   // `ignore_pp` accepts an arbitrary set of arguments, but simply ignores them and returns 0.
   // This can be used in conjunction with parameter pack expansions to execute actions
   // on the packs, by phrasing them as argument expressions. For example, `CxVecK.h`
   // uses a construction like this to implement arithmetic assignment ops for vectors:
   //
   // ```
   //    template<class T, size_t N, cx::index_t... Is>
   //    constexpr auto operator *= (TVecK<T,N,cx::index_list_t<Is...> > &lhs, T scale) -> TVecK<Ta,N,cx::index_list_t<Is...> >& {
   //        return (cx::ignore_pp((lhs[Is] *= scale, 0)...), lhs);
   //    }
   //```
   // WARNING: when using this pattern, the argument sub-expressions may be
   // evaluated in any order the compiler sees fit! C++ makes no guarantee
   // regarding the relative sequencing of function argument evaluations!
   template<class... Ts> int constexpr ignore_pp(Ts...) { return 0; }


} _Namespace_CxMeta_End


_Namespace_CxMeta_Begin {

   template<class... T> class MagicTypeRevealer;
   #ifdef INCLUDE_ABANDONED
   // template<class T> class MagicTypeRevealer;
   //
   // int bla() {
   //    auto Is = MakeIndexList<5>();
   //    MagicTypeRevealer<decltype(Is)>();
   // }
   #endif // INCLUDE_ABANDONED


   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: C++11 abbreviations for type transforms
   // ──────────────────────────────────────────────────────────────────────────


   // try to get the actual element type hidden behind something returns as reference,
   // by stripping off reference-ness and const-ness. Of course, if the actual type
   // indeed *is* const, this will remove it, too.
   // (note: there is a std::remove_cvref, but only in >= C++20; I call this a bit
   // differently to reduce the chances of ambiguity with the std version if someone
   // DOES use >= C++20. I did not think about this very deeply, but atm it seems
   // quite possible that multiple variants could get into competition due to ADL
   // or using namespace ...;
   template<class FDeclType>
   using without_cvref_t = typename std::remove_cv<typename std::remove_reference<FDeclType>::type>::type;

   // no std::remove_reference_t in C++11...
   template<class FDeclType>
   using without_ref_t = typename std::remove_reference<FDeclType>::type;

   // attempts to identify the type of the actual objects referenced by iterator
   // range [sequence.begin(), sequence.end()]
   template<class FSequence>
   using sequence_element_type_t = without_cvref_t<decltype(*std::declval<FSequence>().begin())>;


   // a <= C++11 shortcut for the latter std library's std::enable_if_t
   template<bool B, class T = void>
   using enable_if_t = typename std::enable_if<B,T>::type;

   // as in C++14...
   template<bool B, class T, class F>
   using conditional_t = typename std::conditional<B,T,F>::type;

   // helpers for SFINAE-based template specialization
   template<bool B>     struct void_if {};
   template<>           struct void_if<true> { using type = void; };

   template<bool... Bn> struct void_if_any : public void_if<any_pp(Bn...)>::type {};
   template<bool... Bn> struct void_if_all : public void_if<all_pp(Bn...)>::type {};
   template<bool B>     using void_if_t = typename void_if<B>::type;
   template<bool B>     using void_if_not_t = typename void_if<!B>::type;
   template<bool... Bn> using void_if_any_t = typename void_if<any_pp(Bn...)>::type;
   template<bool... Bn> using void_if_all_t = typename void_if<all_pp(Bn...)>::type;

   template<class T, bool B>     struct type_if {};
   template<class T>             struct type_if<T,true> { using type = T; };

   template<class T, bool... Bn> struct type_if_any : public type_if<T, any_pp(Bn...)>::type {};
   template<class T, bool... Bn> struct type_if_all : public type_if<T, all_pp(Bn...)>::type {};
   template<class T, bool B>     using type_if_t = typename type_if<T, B>::type;
   template<class T, bool B>     using type_if_not_t = typename type_if<T, !B>::type;
   template<class T, bool... Bn> using type_if_any_t = typename type_if<T, any_pp(Bn...)>::type;
   template<class T, bool... Bn> using type_if_all_t = typename type_if<T, all_pp(Bn...)>::type;

#ifdef INCLUDE_ABANDONED
//    template<bool... Bn> using void_if_any_t = typename void_if_any<Bn...>::type;
//    template<bool... Bn> using void_if_all_t = typename void_if_all::type;

/*
//    template<> struct void_if_any<> {};
//    template<bool... _Bn> struct void_if_any<true, _Bn...> { using type = void; }
//    template<bool... _Bn> struct void_if_any<false,_Bn...> : public typename void_if_any<_Bn...>::type {}
//    template<typename _B1, typename... _Bn> struct void_if_any<_B1, _Bn...> : public std::conditional<_B1::value, _B1, OrT<_B1, _Bn...> >::type {};

   /// `or_t<C1, C2,...>` evaluates to `void` if any of the integral bool constants C1,C2,...
   /// are `true`, and to SFINAE otherwise.
   template<bool... _Conds>
   using OrC = typename OrT<typename enable_if<_Conds,FResult>...>::type;

   template<typename...> struct AndT;
   template<> struct AndT<> : public std::true_type {};
   template<typename _B1, typename... _Bn> struct AndT<_B1, _Bn...> : public std::conditional<_B1::value, _B1, AndT<_B1, _Bn...> >::type {};

   template<typename B>
   struct NotT : public std::integral_contant<bool, !B::value> {};

   template<typename... Conds>
   using _require = typename enable_if<__and_<_Cond...>::value>::type;*/
#endif // INCLUDE_ABANDONED


#ifdef INCLUDE_ABANDONED
//    // unwrap_decay_t recycled from:
//    //
//    //     https://en.cppreference.com/w/cpp/utility/tuple/make_tuple
//    //
//    // It just runs through the normal std::decay, except if doing so results in
//    // a reference wrapper. In that case, the reference is unwrapped and an
//    // actual ref is created.
//    //
//    // This is process is what make_tuple() uses to decide on output types. And
//    // that is probably rather important when it comes to matters like not
//    // creating dangling r-value references...
//    template <class T> struct unwrap_ref{ using type = T; };
//    template <class T> struct unwrap_ref<std::reference_wrapper<T> > { using type = T&; };
//    template <class T> using unwrap_ref_t = typename unwrap_ref<T>::type;
//
//    template <class T> using decay_and_unwrap_refs_t = unwrap_ref_t<decay_t<T> >;
//    // ^-- >= C++20 has std::unwrap_ref_decay_t
//
// ^-- moved to CxMeta.h
#endif // INCLUDE_ABANDONED

   template <class T> using decay_t = typename std::decay<T>::type;

   namespace iter {
      namespace detail {
         template<class T>
         struct _TGetSizeHelper {
            static auto constexpr size(T const &seq) -> decltype(seq.size()) { return seq.size(); }
         };
         template<class T, std::size_t N>
         struct _TGetSizeHelper<T[N]> {
            static std::size_t constexpr size(T const *) { return N; }
         };
      }

      /// meant as C++11-level substitute for std::size.
      template<class T> std::size_t constexpr _SequenceSizeQ(T const &seq) noexcept {
         return detail::_TGetSizeHelper<cx::without_cvref_t<T> >::size(seq);
      }

      using std::begin; using std::end;
      // ^-- there is no other way of doing that, is there?
      template<class FSeq> auto _SequenceBeginQ(FSeq &seq) noexcept -> cx::without_ref_t<decltype(begin(seq))> { return begin(seq); }
      template<class FSeq> auto _SequenceEndQ(FSeq &seq) noexcept -> cx::without_ref_t<decltype(end(seq))> { return end(seq); }
      template<class FSeq> using _SequenceIterT = decltype(begin(std::declval<FSeq>()));
      template<class FSeq> using _SequenceDerefIterT = decltype(*begin(std::declval<FSeq>()));
      // ^-- note: for most types of sequences (in particular, actual containers)
      //     this will result in a reference type!
      //     ...this is why I did not call that ValueT.
      template<class FSeq> using _SequenceValueT = typename std::remove_reference<decltype(*begin(std::declval<FSeq>()))>::type;
   } // namespace iter

#ifdef INCLUDE_ABANDONED
//    template<class T, std::size_t N>
//    std::size_t constexpr sequence_size(const T (&array)[N]) noexcept { return N; }
//
//    template<class C> constexpr auto sequence_size(const C& c) -> decltype(c.size()) { return c.size(); }
#endif // INCLUDE_ABANDONED


   // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   // ▒▒ NEXT: C++11 abbreviations for tuple-interface sizes and types/values
   // ──────────────────────────────────────────────────────────────────────────

//    template<class FTupleLike>
//    inline constexpr size_t tuple_size_v = std::tuple_size<FTupleLike>::value;

#if __cplusplus >= 201703L
   template<class FTupleLike> inline constexpr size_t _NumElemQ = std::tuple_size<FTupleLike>::value;
   template<class FTupleLike> inline constexpr size_t _ElemNumQ = std::tuple_size<FTupleLike>::value;
   // ^-- can't decide. Also count vs num...?
#endif // __cplusplus >= 201703L
   template<class FTupleLike> using _NumElemT = typename std::tuple_size<without_cvref_t<FTupleLike> >;
   template<class FTupleLike> using _ElemNumT = typename std::tuple_size<without_cvref_t<FTupleLike> >;

   template<size_t iElem, class FTupleLike> using _ElemTypeQ = typename std::tuple_element<iElem, FTupleLike>::type;

//    template<size_t iElem, class FTupleLike>
//    auto _ElemQ(FTupleLike &&obj) -> decltype(std::get<iElem>(std::forward(obj))) { return std::get<iElem>(std::forward(obj)); };
   // ^-- _ElemValueQ?


//    template<size_t iElem, class FTupleLike>
//    auto _ElemFwdQ(FTupleLike &&obj) -> decltype(get<iElem>(std::forward<FTupleLike>(obj))) { return get<iElem>(std::forward<FTupleLike>(obj)); }
   template<size_t iElem, class FTupleLike>
   auto _ElemFwdQ(FTupleLike &&obj) -> _ElemTypeQ<iElem, FTupleLike>&& { using std::get; return get<iElem>(std::forward<FTupleLike>(obj)); }
   // ^--- I am still not good at modern C++... and I wonder if there are ever
   // cases where I should *not* use get<...>(forward)? needs some more
   // thinking, I guess...
   // maybe cross check: https://en.cppreference.com/w/cpp/utility/forward
   template<size_t iElem, class FTupleLike>
   auto _ElemValQ(FTupleLike const &obj) -> _ElemTypeQ<iElem, FTupleLike> { using std::get; return get<iElem>(obj); }
//    auto _ElemValQ(FTupleLike const &obj) -> decltype(get<iElem>(obj)) { using std::get; return get<iElem>(obj); }

   template<class FTupleLike, class FIndexList = index_list_of_size_t<_NumElemT<FTupleLike>::value > >
   auto _MakeElemIndexList(FTupleLike const &) -> FIndexList { return FIndexList(); }


#ifdef INCLUDE_ABANDONED
/*
   // FIXME: I think the full implementations of all this stuff are in ~/dev/cx/tests/list_indices, or somewhere around.

   namespace iter_detail {
      using std::get;
      template<class FTupleLike, size_t... Is>
      auto _AsTupleOfRefsImpl(FTupleLike &&obj, cx::index_list_t<Is...>)
      -> decltype(std::forward_as_tuple(get<Is>(std::forward<FTupleLike>(obj))...))
      {
      //    return {FixedSizeArray[Is]...};
         return std::forward_as_tuple(get<Is>(std::forward<FTupleLike>(obj))...);
      }

      template<class FTupleLike, size_t... Is>
      auto _AsTupleOfValsImpl(FTupleLike &&obj, cx::index_list_t<Is...>)
      -> std::tuple<std::decay<decltype(get<Is>(obj))>...>
      {
         return {get<Is>(obj)...};
      }
   }

   // make an std::tuple of forwarding references for the elements of obj (e.g., a fixed size array)
   template<class FTupleLike>
   auto _AsTupleOfRefs(FTupleLike &&obj)
   -> decltype(iter_detail::_AsTupleOfRefsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj)))
   {
      return iter_detail::_AsTupleOfRefsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj));
   }

   // make an std::tuple of actual values for the input elements (i.e., not reference or
   // forwarding references)
   template<class FTupleLike>
   auto _AsTupleOfVals(FTupleLike &&obj)
   -> decltype(iter_detail::_AsTupleOfValsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj)))
   {
      return iter_detail::_AsTupleOfValsImpl<FTupleLike>(std::forward<FTupleLike>(obj), _MakeElemIndexList(obj));
   }

   */
#endif // INCLUDE_ABANDONED

#ifdef INCLUDE_ABANDONED
//    // see here: for more clever versions (for types):
//    //
//    //  https://ldionne.com/2015/11/29/efficient-parameter-pack-indexing/
//    //
//    // And the comment that it apparently could be implemented
//    // in a few hours into an existing compiler...  there apparently are even clang and g++ compiler extensions!
//    namespace detail {
//       using std::get;
//       template<size_t N, class FTupleLike, index_t... Is>
//       auto forward_first_n_impl(FTupleLike &&objs, index_list_t<Is...>)
//       -> decltype(std::forward_as_tuple(get<Is>(objs)...))
//       {
//          return std::forward_as_tuple(get<Is>(objs)...);
//       }
//    }
//
//    template<size_t N, class FTupleLike>
//    auto forward_first_n(FTupleLike &&objs)
//    -> decltype(detail::forward_first_n_impl<N>(std::forward<FTupleLike>(objs), index_list_of_size_t<N>()))
//    {
//       static_assert(N <= std::tuple_size<FTupleLike>::value, "forward_first_n<N>(T obj): type T does not have N elements to forward.");
//       return detail::forward_first_n_impl<N>(std::forward<FTupleLike>(objs), index_list_of_size_t<N>());
//    }
// ^-- moved to CxMeta.h
#endif // INCLUDE_ABANDONED

} _Namespace_CxMeta_End


_Namespace_CxMeta_Begin {
   /// auxiliary class with a member enum {value} which evaluates to `true` if all types
   /// in Args... are equal to RefType, and false otherwise.
   template<class RefType, class... Args>
   struct all_types_equal;

   template<class RefType>
   struct all_types_equal<RefType> {
      enum { value = true };
   };

   template<class RefType, class Arg0, class... Args>
   struct all_types_equal<RefType, Arg0, Args...> {
      enum { value = std::is_same<RefType, Arg0>::value && all_types_equal<RefType, Args...>::value };
   };


   template<class... Args>
   using all_types_equal_v = typename all_types_equal<Args...>::value;

#ifdef INCLUDE_ABANDONED
//    namespace detail {
//       // typedefs ::type, which is either std::array<Arg,N> if all argument data
//       // types are identical, or std::tuple<Args...> if not.
//       template<class... Args> struct TPacket;
//
//       // note: gets instanciated with cx::without_cvref_t<Args>..., not Args... directly!
//       // (see _SelectAxisLayoutType typedef below)
//       template<class Arg0, class... Args>
//       struct TPacket<Arg0, Args...> {
//          using type = cx::conditional_t<
//             cx::all_types_equal<Arg0, Args...>::value,
//                std::array<Arg0, 1 + sizeof...(Args)>,
//                std::tuple<Args...>
//             >;
//       };
//       // add a special case for zero-argument packet. That needs to be a
//       // tuple, because zero-length arrays are not allowed.
//       template<> struct TPacket<> {
//          using type = std::tuple<>;
//       };
//    }
//
//    /// Resolves to either std::tuple<Args...> or std::array<Arg,N>, depending on whether all
//    /// class arguments are equal. The types are taken literally---i.e., if the entries in Args...
//    /// have any reference or cv qualifiers, so will the entries in the returned tuple.
//    template<class... Args> using decltype_packet_t = typename detail::TPacket<Args...>::type;
//
//    /// Resolves to either a std::tuple or a std::array suitable for storing copies of the values
//    /// of the expressions in Args...; that is, if the values have any reference or cv qualifiers,
//    /// those are removed.
//    template<class... Args> using value_packet_t = typename detail::TPacket<cx::without_cvref_t<Args>...>::type;
//
//    template<class... Args> auto make_decltype_packet(Args... args) -> decltype_packet_t<Args...> { return {args...}; }
//    template<class... Args> auto make_value_packet(Args... args) -> value_packet_t<Args...> { return {args...}; }
// ^-- moved to CxMeta.h
#endif // INCLUDE_ABANDONED
} _Namespace_CxMeta_End



/*
/// Defines a minimal iterator interface to allow use in range-based for loops
/// and generic container printing mechanism. This one is for objects which have
/// a working [] operator to access elements, but not necessarily much else.
/// Examples are raw-pointer based arrays (with non-unit stride, otherwise one
/// can just use the pointers directly).
/// Usage: Add:
///
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,const_iterator,const)
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,iterator,)
///
/// to class declaration to make both const and non-const iterators which wrap
/// calls to operator [].
#define CX_IMPLEMENT_ITERATOR_IndexInBracket(HostType,IteratorType,ConstOrNothing) \
   struct IteratorType { \
      HostType ConstOrNothing \
         *m_pArrayObj; \
      size_t \
         m_Index; \
      IteratorType &operator ++() { ++m_Index; return *this; } \
      IteratorType &operator --() { --m_Index; return *this; } \
      bool operator != (IteratorType const& other) const { \
         return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      bool operator == (IteratorType const& other) const { \
         return this->m_pArrayObj == other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      auto operator *() -> decltype((*m_pArrayObj)[m_Index]) { return (*m_pArrayObj)[m_Index]; } \
   }; \
   IteratorType begin() ConstOrNothing { return IteratorType{this, 0}; } \
   IteratorType end() ConstOrNothing { return IteratorType{this, size()}; }
   // ^--- note: == operator not strictly required.
//       typename std::decay<decltype(std::declval<HostType>()[0])>::type
//       ConstOrNothing &operator *() { return (*m_pArrayObj)[m_Index]; }
^-- moved to CxDefs.h for the moment.
*/

#ifdef INCLUDE_ABANDONED
   // RANT: ...seriously what were the C++ standards guys thinking when they designed this??!
   //
   // - No way to access parameter packs by index, neither types nor values
   //   (it can be done by meta-programming hacks like here, but do we really want that?)
   //
   //   - No way to create a correspondence to parameter packs inside a function.
   //     They indeed can *only* be made or used with actual parameters.
   //     So just about *every single function* we might be dealing with needs
   //     to be split into one function to generate an index list for a parameter pack,
   //     and another functions which receives it as argument, for the simple reason
   //     that there is no other way to access them.
   //
   //   - In C++17 there are fold expressions... but only for parameter packs (e.g.,
   //     not any of the other things implementing get<...>), and only based on
   //     a very small number of predefined arithmetic operators (not generic unary
   //     or binary functions)... WHY?
   //
   // - The way tuples work...
   //   - instead of specifying a straight-forward compiler-level implementation
   //     (which they needed ANYWAY, for the parameter packs! Because tuples
   //     are exactly what those are), the entire thing is meant to be based on some
   //     insanity template meta-programmng mess.
   //
   //   - And the syntax is just stupid. Admittedly, my >= C++11 experience is still
   //     limited, but in my trials so far, these things have been mostly cumbersome
   //     (not the least because of all the forwarding/rvalue ref stuff).
   //     In most cases simply *not* using tuples had been the better option!
   //
   //   - e.g., why is there no operator []? There are const-based overloads,
   //     why not constexpr based ones, too? instead, it's creating an own local
   //     namespace, `using std::get` inside this namespace (because you can't
   //     *just* use std::get if you want to be compatible, because other `get`
   //     methods need to be resolved via ADL), then `get<...>`(t) to access
   //     elements). Was there really *no way at all* to add an operator []?
   //     Yes, different tuple elements may have different types, but so have
   //     the expressions get<n>(T). I cannot believe that this couldn't easily
   //     be fixed with a bit of compiler magic (and tuples really *SHOULD*
   //     be implemented in the compiler).
   //
   //   - Tuples are super easy to use unsafely, including when using them for
   //     their intended core usage scenarios, and when following standard
   //     recommended usage patterns to the letter. I already managed to create
   //     several situations in trial programs in which even with all compiler
   //     warnings and memory sanitizing etc enabled, I clearly could access
   //     variables beyond their end of lifetime (admittedly, I explicitly
   //     *made* the programs to check for situations in which this would
   //     arise).
   //
   // - basic aggregates can be used in aggregate initialization ({}) and
   //   structured binding, including simple structs like `{int a; int b;}`.
   //   Very nice---in principle. But to make this practically usable, it is
   //   clear that the designers of the language must have recognized that in
   //   many situations:
   //
   //   - structs have a well-defined number and order of arguments
   //
   //   - and the users of the language can recognize their meaning even without
   //     element names (just like they can do so for function arguments...)
   //
   //   So far so good. So... why then is there no way to access these elements in
   //   a generic way by index in *any other* situation than initialization or
   //   structured binding?! Why is there no std::get<n> for aggregates?
#endif // INCLUDE_ABANDONED


#endif // CX_META_LIGHT_H
