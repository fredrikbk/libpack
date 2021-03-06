
define void @packer(i8* %inbuf, i32 %count, i8* %outbuf) nounwind {
entry:
  %out = ptrtoint i8* %outbuf to i64
  %in = ptrtoint i8* %inbuf to i64
  br label %outerloop

outerloop:                                        ; preds = %afterinner, %entry
  %out1 = phi i64 [ %out, %entry ], [ %0, %afterinner ]
  %in1 = phi i64 [ %in, %entry ], [ %nextin1, %afterinner ]
  %i = phi i32 [ 0, %entry ], [ %nexti, %afterinner ]
  %0 = add i64 %out1, 1280
  br label %innerloop

innerloop:                                        ; preds = %innerloop, %outerloop
  %out2.1 = phi i64 [ %out1, %outerloop ], [ %out2.5, %innerloop ]
  %in2.1 = phi i64 [ %in1, %outerloop ], [ %in2.5, %innerloop ]

  %out2_addr.1 = inttoptr i64 %out2.1 to <40 x i8>*
  %in2_addr.1 = inttoptr i64 %in2.1 to <40 x i8>*
  %out2.2 = add i64 %out2.1, 40
  %in2.2 = add i64 %in2.1, 1440

  %out2_addr.2 = inttoptr i64 %out2.2 to <40 x i8>*
  %in2_addr.2 = inttoptr i64 %in2.2 to <40 x i8>*
  %out2.3 = add i64 %out2.2, 40
  %in2.3 = add i64 %in2.2, 1440

  %out2_addr.3 = inttoptr i64 %out2.3 to <40 x i8>*
  %in2_addr.3 = inttoptr i64 %in2.3 to <40 x i8>*
  %out2.4 = add i64 %out2.3, 40
  %in2.4 = add i64 %in2.3, 1440

  %out2_addr.4 = inttoptr i64 %out2.4 to <40 x i8>*
  %in2_addr.4 = inttoptr i64 %in2.4 to <40 x i8>*
  %out2.5 = add i64 %out2.4, 40
  %in2.5 = add i64 %in2.4, 1440

  %in2_addr.2.prefetch = inttoptr i64 %in2.2 to i8*
  call void @llvm.prefetch(i8* %in2_addr.2.prefetch, i32 0, i32 0, i32 1)
  %bytes.1 = load <40 x i8>* %in2_addr.1, align 1
  store <40 x i8> %bytes.1, <40 x i8>* %out2_addr.1, align 1

  %in2_addr.3.prefetch = inttoptr i64 %in2.3 to i8*
  call void @llvm.prefetch(i8* %in2_addr.3.prefetch, i32 0, i32 0, i32 1)
  %bytes.2 = load <40 x i8>* %in2_addr.2, align 1
  store <40 x i8> %bytes.2, <40 x i8>* %out2_addr.2, align 1

  %in2_addr.4.prefetch = inttoptr i64 %in2.4 to i8*
  call void @llvm.prefetch(i8* %in2_addr.4.prefetch, i32 0, i32 0, i32 1)
  %bytes.3 = load <40 x i8>* %in2_addr.3, align 1
  store <40 x i8> %bytes.3, <40 x i8>* %out2_addr.3, align 1

  %bytes.4 = load <40 x i8>* %in2_addr.4, align 1
  store <40 x i8> %bytes.4, <40 x i8>* %out2_addr.4, align 1

  %innercond = icmp eq i64 %out2.5, %0
  br i1 %innercond, label %afterinner, label %innerloop

afterinner:                                       ; preds = %innerloop
  %nextin1 = add i64 %in1, 44680
  %nexti = add i32 %i, 1
  %outercond = icmp eq i32 %nexti, %count
  br i1 %outercond, label %afterouter, label %outerloop

afterouter:                                       ; preds = %afterinner
  ret void
}

declare void @llvm.prefetch(i8* nocapture, i32, i32, i32) nounwind
