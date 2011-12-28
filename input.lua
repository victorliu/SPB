S = SPB.NewBandSolver{
	Lattice = {
		{1,0},
		{0,1}
	},
	Polarization = 'E'
}
scale = 2*math.pi


S:SetOptions{
	NumBands = 30,
	Resolution = {10, 10},
	TargetFrequency = 0. * scale,
	Tolerance = 1e-7,
	Verbosity = 0
}

S:AddMaterial{
	Name = 'Si',
	EpsilonInf = {12,0}
}

S:AddMaterial{
	Name = 'Diel',
	EpsilonInf = {2.3646,0}
}
S:AddMaterial{
	Name = 'Silver',
	EpsilonInf = {2.3646,0}
}
S:AddMaterialLorentzPole{
	Material = 'Silver',
	Omega0 = 0.4593 * scale,
	Gamma = 0,
	OmegaP = 0.1676 * scale
}
S:AddMaterialLorentzPole{
	Material = 'Silver',
	Omega0 = 0.5434 * scale,
	Gamma = 0,
	OmegaP = 0.3293 * scale
}
S:AddMaterialLorentzPole{
	Material = 'Silver',
	Omega0 = 0,
	Gamma = 0,
	OmegaP = 0.6253 * scale
}


S:SetRectangle{
	Material = 'Silver',
	Center = {0,0},
	Halfwidths = {0.225,0.225},
	Angle = 0
}
--[[
S:OutputEpsilon{
	Resolution = {5,5},
	Filename = "epsilon.txt",
	Format = 'gnuplot'
}
os.exit(0)
]]
for kx = 0.01,0.5,0.01 do
	S:SolveK{kx,0}

	io.stdout:write(kx)
	for i,freq in ipairs(S:GetFrequencies()) do
		if freq[1] >= 0 then
			io.stdout:write('\t' .. freq[1]/scale)
		end
	end
	io.stdout:write('\n')
	io.stdout:flush()
end


--[[
band1 = S:GetBand(1);
band2 = S:GetBand(2);
dot = S:InnerProduct(band1,band2)
band1:OutputHDF{
	BaseName = 'field_'
}
]]
