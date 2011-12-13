S = SPB.NewBandSolver{
	Lattice = {
		{1,0},
		{0,1}
	},
	Polarization = 'E'
}

S:SetOptions{
	NumBands = 10,
	Resolution = {10,10},
	TargetFrequency = 0.1,
	Tolerance = 1e-7,
	Verbosity = 0
}

S:AddMaterial{
	Name = 'Si',
	EpsilonInf = {12,0}
}

S:AddMaterial{
	Name = 'Metal',
	EpsilonInf = {1,0}
}
S:AddMaterialLorentzPole{
	Material = 'Metal',
	Omega0 = 0.5,
	Gamma = 0,
	OmegaP = 1
}


S:SetRectangle{
	Material = 'Metal',
	Center = {0,0},
	Halfwidths = {0.25,0.25},
	Angle = 0
}
--[[
S:OutputEpsilon{
	Resolution = {32,32},
	Filename = "epsilon.txt",
	Format = 'gnuplot'
}
]]
for kx = 0,0.5,0.01 do
	S:SolveK{kx,0}

	io.stdout:write(kx)
	for i,freq in ipairs(S:GetFrequencies()) do
		if freq[1] >= 0 then
			io.stdout:write('\t' .. freq[1])
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
