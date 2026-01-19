import csv
import matplotlib.pyplot as plt

TOP_N = 5
RXMIN_DBM = -25.0 

def read_ont_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            row["ont_id"] = int(row["ont_id"])
            row["total_dist_km"] = float(row["total_dist_km"])
            row["total_loss_db"] = float(row["total_loss_db"])
            row["rx_dbm"] = float(row["rx_dbm"])
            row["margin_db"] = float(row["margin_db"])
            rows.append(row)
    return rows

def read_splitter_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            row["ratio"] = int(row["ratio"])
            row["ont_count"] = int(row["ont_count"])
            row["ok_count"] = int(row["ok_count"])
            row["fail_count"] = int(row["fail_count"])
            row["down_count"] = int(row["down_count"])
            row["avg_rx_dbm"] = float(row["avg_rx_dbm"])
            row["avg_loss_db"] = float(row["avg_loss_db"])
            row["worst_rx_dbm"] = float(row["worst_rx_dbm"])
            rows.append(row)
    return rows

def main():
    onts = read_ont_csv("ont_results.csv")
    splits = read_splitter_csv("splitter_results.csv")

    # ---- priprema figure: 2 reda x 2 stupca ----
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.suptitle("FTTH/GPON analiza optičke mreže", fontsize=14, fontweight="bold")

    # 1) Histogram RX snage (bez DOWN)
    rx_vals = [o["rx_dbm"] for o in onts if o["status"] != "DOWN"]
    axes[0, 0].hist(rx_vals, bins=12, alpha=0.75)
    axes[0, 0].set_title("Distribucija RX snage (bez DOWN)")
    axes[0, 0].set_xlabel("RX snaga [dBm]")
    axes[0, 0].set_ylabel("Broj ONT-ova")
    axes[0, 0].grid(True, linestyle="--", alpha=0.6)

    # 2) RX snaga vs udaljenost (boje po statusu) + RXmin linija
    dist_ok   = [o["total_dist_km"] for o in onts if o["status"] == "OK"]
    rx_ok     = [o["rx_dbm"] for o in onts if o["status"] == "OK"]

    dist_fail = [o["total_dist_km"] for o in onts if o["status"] == "FAIL"]
    rx_fail   = [o["rx_dbm"] for o in onts if o["status"] == "FAIL"]

    dist_down = [o["total_dist_km"] for o in onts if o["status"] == "DOWN"]
    rx_down   = [o["rx_dbm"] for o in onts if o["status"] == "DOWN"]

    axes[0, 1].scatter(dist_ok, rx_ok, alpha=0.7, label="OK")
    axes[0, 1].scatter(dist_fail, rx_fail, alpha=0.7, label="FAIL")
    axes[0, 1].scatter(dist_down, rx_down, alpha=0.7, label="DOWN")

    # RXmin horizontalna linija
    axes[0, 1].axhline(RXMIN_DBM, linestyle="--", linewidth=2, label=f"RXmin = {RXMIN_DBM:.1f} dBm")

    axes[0, 1].set_title("RX snaga u odnosu na udaljenost (po statusu)")
    axes[0, 1].set_xlabel("Ukupna udaljenost [km]")
    axes[0, 1].set_ylabel("RX snaga [dBm]")
    axes[0, 1].grid(True, linestyle="--", alpha=0.6)
    axes[0, 1].legend()

    # 3) FAIL + DOWN po splitteru
    names = [s["name"] + f" (1:{s['ratio']})" for s in splits]
    bad = [s["fail_count"] + s["down_count"] for s in splits]

    axes[1, 0].bar(names, bad, alpha=0.8)
    axes[1, 0].set_title("Neispravni ONT po splitteru")
    axes[1, 0].set_xlabel("Splitter")
    axes[1, 0].set_ylabel("FAIL + DOWN")

    axes[1, 0].tick_params(axis="x", rotation=45)
    for lbl in axes[1, 0].get_xticklabels():
        lbl.set_ha("right")

    axes[1, 0].grid(True, axis="y", linestyle="--", alpha=0.6)


    # 4) TOP N najgorih ONT-ova po margini
    onts_sorted = sorted(onts, key=lambda o: o["margin_db"])
    worst = onts_sorted[:TOP_N]

    worst_labels = [f"ONT {o['ont_id']}" for o in worst]
    worst_margins = [o["margin_db"] for o in worst]
    worst_rx = [o["rx_dbm"] for o in worst]

    axes[1, 1].bar(worst_labels, worst_margins, alpha=0.85)
    axes[1, 1].set_title(f"TOP {TOP_N} najgorih ONT-ova (po margini)")
    axes[1, 1].set_xlabel("ONT")
    axes[1, 1].set_ylabel("Margina [dB] (RX - RXmin)")
    axes[1, 1].grid(True, axis="y", linestyle="--", alpha=0.6)

    # tekst: pošto su margine negativne, stavi labelu malo IZNAD dna stupca (prema gore)
    for i, (m, rx) in enumerate(zip(worst_margins, worst_rx)):
        axes[1, 1].text(
            i, m + 1.0,              # +1 dB prema gore (unutra u stupac)
            f"RX {rx:.1f}",
            ha="center",
            va="bottom",
            fontsize=9,
            clip_on=True
        )
    
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.12, top=0.90, wspace=0.25, hspace=0.35)
    plt.show()

if __name__ == "__main__":
    main()
