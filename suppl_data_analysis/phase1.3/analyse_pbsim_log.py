import statistics

from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader

if __name__ == '__main__':
    retd = {"I": [], "D": [], "M": [], "S": []}
    with get_tqdm_line_reader("pbsim.log") as reader:
        for line in reader:
            if line.startswith("substitution rate"):
                retd["S"].append(float(line.split(":")[1].strip()))
            elif line.startswith("insertion rate"):
                retd["I"].append(float(line.split(":")[1].strip()))
            elif line.startswith("deletion rate"):
                retd["D"].append(float(line.split(":")[1].strip()))
            elif line.startswith("read accuracy mean"):
                retd["M"].append(float(line.split(":")[1].split("(")[0].strip()))
        print("Events:", {k: str(round(statistics.mean(v) * 100, 2)) + "%" for k, v in retd.items()})
