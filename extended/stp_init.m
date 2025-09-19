function STP = stp_init(Nedges, Reg)
STP.R = ones(Nedges,1);
STP.u = Reg.STP.U0*ones(Nedges,1);
STP.D = Reg.STP.D; STP.F = Reg.STP.F;
end
