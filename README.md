# FTTH / GPON Network Simulator

C · Data Structures Project

# 📌 About

FTTH/GPON network simulator developed as part of the Data Structures course.
The project simulates an optical access network and analyzes signal quality for each user (ONT).

🧠 Data Structures

🌳 N-ary tree (child / sibling) – models OLT → Splitter → ONT hierarchy

📊 Aggregation struct – collects statistics during recursion

📋 Dynamic array – stores per-splitter results

# 🚀 Features

Reads network topology from a text file

Recursive network traversal

Calculates optical losses and RX power

Determines ONT status (OK / FAIL / DOWN)

Aggregates statistics per splitter

Exports results to CSV and report files

Identifies TOP N worst ONTs by margin

# 📈 Visualization

Python script generates graphs:

RX power distribution

RX vs distance (with RXmin line)

FAIL + DOWN per splitter

TOP N worst ONTs

# 🎯 Focus

Language: C

Course: Data Structures

Concepts: Trees, Recursion, Dynamic Memory, Algorithms
