#!/usr/bin/env python3
from table import TableResults
from table import SortList, SortListEntry
from exceptions import DuplicateKeyError

tr = TableResults("software X", "6.6.6", "20200128", "runsoft -i arg",
                  "sampleX")

slentry = SortListEntry(1, 2)
print("Entry: {}".format(slentry))

list1 = [1, 2, 3, 4]
list2 = ["d", "c", "b", "a"]
list3 = [1, 2]

slist0 = SortList()
slist0.append("k1", "v1")
slist0.append("k2", "v2")
print("slist0:\n{}".format(slist0))

slist = SortList(unique_list=list1, val_list=list2)

print(slist)

try:
    SortList(unique_list=list1, val_list=list3)
except IndexError as e:
    print("Caught expected IndexError:\n{}".format(e))

slist.sort()
print("Sorted list:\n{}".format(slist))

tr.add_table("My table")
mytbl = tr.long["My table"]
mytbl.add_headers(["header1", "header2", "header3"])

mytbl["row2"] = {
    "header1": "a",
    "header2": "b",
    "header3": "c",
}
mytbl["row1"] = {
    "header1": "A",
    "header2": "B",
    "header3": "C",
}
mytbl["row3"] = {
    "header1": 1,
    "header3": 2,
}

print("My table:\n{}".format(tr.long["My table"]))

print("my table sort list:\n{}".format(tr.long["My table"]._sort_list))

print("row headers: {}"
    .format(tr.long["My table"].extract_column("row_header")))
print("header1: {}"
    .format(tr.long["My table"].extract_column("header1")))
print("header2: {}"
    .format(tr.long["My table"].extract_column("header2")))

mytbl["row1"]["header4"] = "X"
print("Added single value with new header:\n{}".format(mytbl))
mytbl["row4"] = {
    "header1": "Y",
    "header3": "Y",
    "header5": "Y",
}
print("Added new row with new header:\n{}".format(mytbl))
mytbl["row1"] = {
    "header1": "AA",
    "header3": "BB",
    "header5": "DD",
}
print("Added data via dict to existing row:\n{}".format(mytbl))
print("\theader1\n{}".format(mytbl.extract_column("header1")))

my_slist = tr.long["My table"].get_sort_list()
print("My slist:\n{}".format(my_slist))
tr.long["My table"].set_sort_key("header1")
my_slist = tr.long["My table"].get_sort_list()
print("My slist (sorted header1):\n{}".format(my_slist))

try:
    my_slist.append(2, 4)
except IndexError as e:
    print("Caught expected IndexError:\n{}".format(e))

print("Print rows in order")
print("\t".join(tr.long["My table"].get_headers()))
for row in tr.long["My table"]:
    row_list = tr.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

tr.add_table("My second table")
mytbl2 = tr.long["My second table"]
mytbl2.add_headers(["header4", "header5", "header6"])
mytbl2.lock_headers = True
mytbl2["row2"] = {
    "header4": "a",
    "header5": 10,
    "header6": "c",
}
mytbl2["row1"] = {
    "header4": "a",
    "header5": 1,
    "header6": "a",
}
mytbl2["row3"] = {
    "header4": "b",
    "header5": 1000,
    "header6": "b",
}
my_slist2 = tr.long["My second table"].get_sort_list()
print("My slist2 (no sort):\n{}".format(my_slist2))
tr.long["My second table"].set_sort_key("header4")
my_slist2 = tr.long["My second table"].get_sort_list()
print("My slist2 (sort on 4):\n{}".format(my_slist2))
tr.long["My second table"].set_sort_key("header6")
my_slist2 = tr.long["My second table"].get_sort_list()
print("My slist2 (sort on 6):\n{}".format(my_slist2))
tr.long["My second table"].set_sort_key("header5")
my_slist2 = tr.long["My second table"].get_sort_list()
print("My slist2 (sort on 5):\n{}".format(my_slist2))

print("Print second table")
print("\t".join(tr.long["My second table"].get_headers()))
for row in tr.long["My second table"]:
    row_list = tr.long["My second table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

print("Tables in tr:")
for name in tr.long:
    print(name)

tr2 = TableResults("software X", "6.6.6", "20200128", "runsoft -i arg",
                   "sampleX")
tr2.add_table(tr.long["My table"])
del(tr.long["My table"])

print("Tables in tr:")
for name in tr.long:
    print(name)
print("Tables in tr2:")
for name in tr2.long:
    print(name)

tr2.merge(tr)

print("-------------------------------")

print("Tables in tr2:")
for name in tr2.long:
    print(name)

tr.add_table("My table")
mytbl = tr.long["My table"]
mytbl.add_headers(["header1", "header2", "header3", "headerA"])
mytbl["row2"] = {
    "header1": "override a",
    "headerA": "new val",
    "header3": "override c",
}
mytbl["row1"] = {
    "header2": "override B"
}
print("Print tr My table:")
for row in tr.long["My table"]:
    row_list = tr.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

print("Tables in tr:")
for name in tr.long:
    print(name)

print("-------------------------------")

del(tr.long["My second table"])

print("Tables in tr2:")
for name in tr2.long:
    print(name)

print("Print tr2 My table (before merge):")
for row in tr2.long["My table"]:
    row_list = tr2.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

try:
    tr2.merge(tr)
except DuplicateKeyError:
    print("Caught expected duplicate error.")

tr2.long["My table"].rename_row("row1", "row10")
tr2.long["My table"].rename_row("row2", "row20")

print("Tables in tr2:")
for name in tr2.long:
    print(name)

print("Print tr2 My table (after rename):")
for row in tr2.long["My table"]:
    row_list = tr2.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

tr2.merge(tr)

print("Print tr2 My table (after merge):")
for row in tr2.long["My table"]:
    row_list = tr2.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))

# tr2.long["My table"].set_sort_key("header1")
print("SortList:\n{}".format(tr2.long["My table"].get_sort_list()))
tr2.long["My table"].set_sort_key("header1")
print("SortList:\n{}".format(tr2.long["My table"].get_sort_list()))
tr2.long["My table"]["row3"]["header1"] = "AA"
print("SortList:\n{}".format(tr2.long["My table"].get_sort_list()))
print("Print tr2 My table (after row3,header1 change):")
for row in tr2.long["My table"]:
    row_list = tr2.long["My table"].get_row_as_list(row, as_txt=True)
    print("\t".join(row_list))
