import ROOT as rt
from shutil import copyfile

in_file_name = 'ewk.root'

operator = 'FT8'

old_file_name = in_file_name
new_file_name = 'ewk_%s.root'%(operator)

copyfile(old_file_name, new_file_name)

yield_ratio_file_name = 'templates/templates_%s_cutoff_none.root'%(operator)
yield_ratio_file = rt.TFile.Open(yield_ratio_file_name, 'READ')
in_file = rt.TFile.Open(new_file_name, 'UPDATE')

selection = 'BLS'
in_histo_name = 'h3_1'

print in_file
for agqc_param in [2.08, 2, 1, 1.94, 0.975, 1.89, 0.5] :
    #in_histo = in_file.Get('%s/%s_%s'%(selection, selection, in_histo_name))
    in_histo = in_file.Get(in_histo_name)

    if in_histo.GetBinContent(in_histo.GetNbinsX() + 2) : print 'WARNING: Input histogram has overflow! Fix this!'

    string_agqc_param = ('%f' % agqc_param).rstrip('0').rstrip('.')
    scaled_histo = in_histo.Clone('%s_%s_rescaled_%s_%s'%(selection, in_histo_name, operator, string_agqc_param))

    for b in range(1, scaled_histo.GetNbinsX()+1) :
        in_yield = scaled_histo.GetBinContent(b)
        func = yield_ratio_file.Get('bin_content_par1_%d'%b)
        print func

        ratio = func.Eval(agqc_param)
        print 'ratio %f'%ratio
        new_yield = in_yield * ratio
        scaled_histo.SetBinContent(b, new_yield)

    scaled_histo.Write()

in_file.Close()
yield_ratio_file.Close() 
