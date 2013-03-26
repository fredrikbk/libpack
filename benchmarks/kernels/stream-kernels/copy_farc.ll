; ModuleID = 'copy_farc.cpp'
target datalayout = "e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v64:64:64-v128:128:128-a0:0:64-s0:64:64-f80:128:128-n8:16:32:64-S128"
target triple = "x86_64-pc-linux-gnu"

%"class.farc::Datatype" = type { i32 (...)**, void (i8*, i32, i8*)*, void (i8*, i32, i8*)*, %"class.llvm::Function"*, %"class.llvm::Function"* }
%"class.llvm::Function" = type opaque
%"class.farc::PrimitiveDatatype" = type { %"class.farc::Datatype", i32, i32, i32 }

@ddt = global %"class.farc::Datatype"* null, align 8
@inbuf = global [8388608 x double] zeroinitializer, align 1
@outbuf = global [8388608 x double] zeroinitializer, align 1

define void @_Z5stagePPcS0_PiPl(i8** nocapture %out, i8** nocapture %in, i32* nocapture %num, i64* nocapture %size) uwtable {
  tail call void @_ZN4farc8DDT_InitEv()
  %1 = tail call noalias i8* @_Znwm(i64 56)
  %2 = bitcast i8* %1 to %"class.farc::PrimitiveDatatype"*
  invoke void @_ZN4farc17PrimitiveDatatypeC1ENS0_13PrimitiveTypeE(%"class.farc::PrimitiveDatatype"* %2, i32 2)
          to label %3 unwind label %5

; <label>:3                                       ; preds = %0
  %4 = bitcast i8* %1 to %"class.farc::Datatype"*
  tail call void @_ZN4farc10DDT_CommitEPNS_8DatatypeE(%"class.farc::Datatype"* %4)
  store i8* bitcast ([8388608 x double]* @inbuf to i8*), i8** %in, align 8, !tbaa !0
  store i8* bitcast ([8388608 x double]* @outbuf to i8*), i8** %out, align 8, !tbaa !0
  store i32 8388608, i32* %num, align 4, !tbaa !3
  store i64 67108864, i64* %size, align 8, !tbaa !4
  ret void

; <label>:5                                       ; preds = %0
  %6 = landingpad { i8*, i32 } personality i8* bitcast (i32 (...)* @__gxx_personality_v0 to i8*)
          cleanup
  tail call void @_ZdlPv(i8* %1) nounwind
  resume { i8*, i32 } %6
}

declare void @_ZN4farc8DDT_InitEv()

declare noalias i8* @_Znwm(i64)

declare void @_ZN4farc17PrimitiveDatatypeC1ENS0_13PrimitiveTypeE(%"class.farc::PrimitiveDatatype"*, i32)

declare i32 @__gxx_personality_v0(...)

declare void @_ZdlPv(i8*) nounwind

declare void @_ZN4farc10DDT_CommitEPNS_8DatatypeE(%"class.farc::Datatype"*)

define void @_Z4copyPvS_i(i8* %out, i8* %in, i32 %num) uwtable {
  %1 = load %"class.farc::Datatype"** @ddt, align 8, !tbaa !0
  tail call void @_ZN4farc8DDT_PackEPvS0_PNS_8DatatypeEi(i8* %in, i8* %out, %"class.farc::Datatype"* %1, i32 %num)
  ret void
}

declare void @_ZN4farc8DDT_PackEPvS0_PNS_8DatatypeEi(i8*, i8*, %"class.farc::Datatype"*, i32)

!0 = metadata !{metadata !"any pointer", metadata !1}
!1 = metadata !{metadata !"omnipotent char", metadata !2}
!2 = metadata !{metadata !"Simple C/C++ TBAA", null}
!3 = metadata !{metadata !"int", metadata !1}
!4 = metadata !{metadata !"long", metadata !1}
