let nint := 1 ;

let {k in int} p1L[k].intervalID := k ;
let {k in int} p1U[k].intervalID := k ;
let {k in int} p2L[k].intervalID := k ;
let {k in int} p2U[k].intervalID := k ;
let {q in int} p1L[q].parameter :=  1;
let {q in int} p1U[q].parameter :=  1;
let {q in int} p2L[q].parameter :=  2;
let {q in int} p2U[q].parameter :=  2;

let p1L[1].branchmode := 0;
let p1L[1].scaling := 0;
let p1L[1].benfit_value := 1;
let p1L[1] := 0.9;
let p1U[1] := 1.1;

let p2L[1] := 0.9;
let p2U[1] := 1;
