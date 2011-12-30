
S = SPB.NewBandSolver{
	Lattice = {
		{1,0},
		{0,1}
	},
	Polarization = 'E'
}

S:SetOptions{
	NumBands = 30,
	Resolution = {10, 10},
	TargetFrequencyRange = {0.2, 0.6},
	ApproximationTolerance = 1e-4,
	Tolerance = 1e-7,
	Verbosity = 0
}

S:AddMaterial{
	Name = 'Si',
	EpsilonInfinity = {12,0}
}

S:AddMaterial{
	Name = 'Silver',
	EpsilonInfinity = 2.3646,
	LorentzPoles = {
		{
			ResonanceFrequency = {0.4593, 0.0587/2},
			PlasmaFrequency = 0.1676
		},
		{
			ResonanceFrequency = {0.5434, 0.115/2},
			PlasmaFrequency = 0.3293
		},
		{
			ResonanceFrequency = {0, 0.0079/2},
			PlasmaFrequency = 0.6253
		}
	}
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
