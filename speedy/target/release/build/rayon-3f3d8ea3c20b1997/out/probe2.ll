; ModuleID = 'probe2.fa4bdb3f8b637a5d-cgu.0'
source_filename = "probe2.fa4bdb3f8b637a5d-cgu.0"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

; probe2::probe
; Function Attrs: nonlazybind uwtable
define void @_ZN6probe25probe17h0ab101fc28b4c1b6E() unnamed_addr #0 {
start:
  ret void
}

attributes #0 = { nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }

!llvm.module.flags = !{!0, !1}
!llvm.ident = !{!2}

!0 = !{i32 8, !"PIC Level", i32 2}
!1 = !{i32 2, !"RtLibUseGOT", i32 1}
!2 = !{!"rustc version 1.85.0 (4d91de4e4 2025-02-17)"}
