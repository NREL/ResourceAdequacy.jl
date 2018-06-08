struct Backcast <: SinglePeriodExtractionMethod end

"""
Backcast extraction:
Returns a `SystemDistribution` with vg/load samples exactly associated with
`dt`. When the values of `systemset.timestamps` are unique this results in a
single vg/load sample included in the returned value.
"""
function extract(params::Backcast, dt::DateTime,
                 systemset::SystemDistributionSet{N1,T1,N2,T2,P}) where {N1,T1,N2,T2,P}

    sample_idxs = findin(systemset.timestamps, [dt])

    return SystemDistribution{N1,T1,P}(systemset.gen_distrs,
                                       systemset.vgsamples[:, sample_idxs],
                                       systemset.interface_labels,
                                       systemset.interface_distrs,
                                       systemset.loadsamples[:, sample_idxs])

end