S = SPB.NewBandSolver{
	Lattice = {
		{1,0},
		{0,1}
	},
	Polarization = 'E'
}
--[[
S:SetOptions{
	NumBands = 10,
	Resolution = 20,
	TargetFrequency = 0.1,
	Tolerance = 1e-7
}

S:AddMaterial{
	Name = 'Metal',
	EpsilonInf = {1,0},
	LorentzPoles = {
		{
			Omega0 = 1,
			Gamma = 0,
			OmegaP = 0
		}
	}
}

S:AddRectangle{
	Material = 'Metal',
	Center = {0,0},
	Halfwidths = {0.25,0.25},
	Angle = 0
}

S:SolveK{0,0}

print(unpack(S:GetFrequencies()));
]]
--[[
band1 = S:GetBand(1);
band2 = S:GetBand(2);
dot = S:InnerProduct(band1,band2)
band1:OutputHDF{
	BaseName = 'field_'
}
]]
