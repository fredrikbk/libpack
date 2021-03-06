
define void @packer_old(i8* %inbuf, i32 %count, i8* %outbuf) nounwind {
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
  %out2 = phi i64 [ %out1, %outerloop ], [ %nextout2, %innerloop ]
  %in2 = phi i64 [ %in1, %outerloop ], [ %nextin2, %innerloop ]
  %out2_addr = inttoptr i64 %out2 to <40 x i8>*
  %in2_addr = inttoptr i64 %in2 to <40 x i8>*
  %bytes = load <40 x i8>* %in2_addr, align 1
  store <40 x i8> %bytes, <40 x i8>* %out2_addr, align 1
  %nextout2 = add i64 %out2, 40
  %nextin2 = add i64 %in2, 1440
  %innercond = icmp eq i64 %nextout2, %0
  br i1 %innercond, label %afterinner, label %innerloop

afterinner:                                       ; preds = %innerloop
  %nextin1 = add i64 %in1, 44680
  %nexti = add i32 %i, 1
  %outercond = icmp eq i32 %nexti, %count
  br i1 %outercond, label %afterouter, label %outerloop

afterouter:                                       ; preds = %afterinner
  ret void
}

