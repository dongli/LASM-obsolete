load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    f1 = addfile("/Users/dongli/LASM/data_for_paper/referee_1_question_1/lasm.deform.case4.Rp3.0.120x60.00000.nc", "r")
    f2 = addfile("/Users/dongli/LASM/data_for_paper/referee_1_question_1/lasm.deform.case4.Rp3.0.120x60.00300.nc", "r")

    q1 = f1->q4/f1->q0
    q2 = f2->q4/f2->q0

    V = f1->vol

    L_1 = sum(abs(q2-q1)*V)/sum(abs(q1)*V)
    L_2 = sqrt(sum((q2-q1)^2*V))/sqrt(sum(q1^2*V))
    L_inf = max(abs(q2-q1))/max(abs(q1))

    print(L_1)
    print(L_2)
    print(L_inf)

end