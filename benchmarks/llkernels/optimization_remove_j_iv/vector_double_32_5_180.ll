
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
  %out2 = phi i64 [ %out1, %outerloop ], [ %nextout2, %innerloop ]
  %in2 = phi i64 [ %in1, %outerloop ], [ %nextin2, %innerloop ]
  %out2_addr = inttoptr i64 %out2 to i8*
  %in2_addr = inttoptr i64 %in2 to i8*
  call void @llvm.memcpy.p0i8.p0i8.i64(i8* %out2_addr, i8* %in2_addr, i64 40, i32 1, i1 false)
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

declare void @llvm.memcpy.p0i8.p0i8.i64(i8* nocapture, i8* nocapture, i64, i32, i1) nounwind
