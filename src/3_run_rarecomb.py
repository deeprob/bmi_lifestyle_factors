import os
import utils as ut

def main(in_dir, in_sub_dirs, combos):

    for isd in in_sub_dirs:
        rinfile = os.path.join(in_dir, isd, "cat_rarecomb_in.csv")
        pentfile = os.path.join(in_dir, isd, "primary.csv")
        sentfile = os.path.join(in_dir, isd, "secondary.csv")
        routfile = os.path.join(in_dir, isd, f"cat_rarecomb_out_{combos}.csv")
        ut.run_rarecomb(rinfile, pentfile, sentfile, routfile, combos)

    rinfile = os.path.join(in_dir, "meta_rarecomb_in.csv")
    pentfile = os.path.join(in_dir, "primary.csv")
    sentfile = os.path.join(in_dir, "secondary.csv")
    routfile = os.path.join(in_dir, f"meta_rarecomb_out_{combos}.csv")
    ut.run_rarecomb(rinfile, pentfile, sentfile, routfile, combos)
    return

if __name__ == "__main__":
    in_dir = "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/"
    in_sub_dirs = ["mental_health", "physical_activity", "diet"]
    combos = [2, 3]
    for combo in combos:
        print(f"Running rarecomb for {combo}")
        main(in_dir, in_sub_dirs, combo)
