using DelimitedFiles
function GenCurveLoad(filename_GenCurve)
    ## source from: https://www.iso-ne.com/isoexpress/web/reports/operations/-/tree/daily-gen-fuel-type
    temp_data = readdlm(filename_GenCurve, ',',header=true)[1][:,4]
    raw_data = zeros(365,24)
    for dd in 1:365
       for tt in 1:24
          if ~isempty(temp_data[24*(dd-1)+tt])
             raw_data[dd,tt] = temp_data[24*(dd-1)+tt]
          end
       end
    end    
    return raw_data
end
