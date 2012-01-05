
S = SPB.NewBandSolver{
	Lattice = {
		{1,0},
		{0,1}
	},
	Polarization = 'E'
}

S:SetOptions{
	Resolution = {20, 20},
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
	Resolution = {20,20},
	Filename = "epsilon.txt",
	Format = 'gnuplot'
}
os.exit(0)
]]

for kx = 0,0.5,0.02 do
	S:SetK{kx,0}
	interval_list = S:GetApproximateFrequencies{
		TargetFrequencyRange = {0.2, 0.6},
		Tolerance = 1e-2
	}

	io.stdout:write(kx)
	for i,interval in ipairs(interval_list) do
		io.stdout:write('\t[' .. interval[1][1] .. ',' .. interval[1][2] .. ']:' .. interval[2]);
	end
	io.stdout:write('\n')
	io.stdout:flush()
	--[[
	freqs,bands = S:GetBandsNear{
		Frequency = 0.32,
		NumBands = 4,
		Tolerance = 1e-7
	}
	pfreqs = S:GetPerturbedFrequencies{
		Frequencies = freqs,
		Bands = bands
	}
	]]
end



--[[
band1 = S:GetBand(1);
band2 = S:GetBand(2);
dot = S:InnerProduct(band1,band2)
band1:OutputHDF{
	BaseName = 'field_'
}
]]
