from lxml import etree
#import argparse
import sys



def read_in_files(flist):
    roots = []
    species = set()
    specie_diff = {}
    rxn_ids = set()
    my_rxn_file = etree.Element("ReactionScheme")

    for fname_xml in flist:
        print(fname_xml)
        tree = etree.parse(fname_xml)
        roots.append(tree.getroot())
        for son in roots[-1]:
            if son.tag == "Specie":
                specie_id = son.attrib["id"]
                specie_kdiff = son.attrib["kdiff"]
                species.add(specie_id)
                if specie_id in specie_diff:
                    if specie_diff[specie_id]!= specie_kdiff:
                        print(specie_id, specie_diff[specie_id], specie_kdiff)
                else:
                    specie_diff[specie_id] = specie_kdiff
                    my_rxn_file.append(son)
            elif son.tag == "Reaction":
                if son.attrib["id"] in rxn_ids:
                    print(son.attrib["id"], " already exists", fname_xml)
                else:
                    rxn_ids.add(son.attrib["id"])
                    my_rxn_file.append(son)
    return my_rxn_file

flist_ER_no_IP3 = [
    "Rxn_module_Ca.xml",
    "Rxn_module_RyR2_CaM.xml",
    "Rxn_module_SERCA2.xml",
    "Rxn_module_CK.xml",
    "Rxn_module_PP2B.xml",
    ]
flist_ER = flist_ER_no_IP3 + ["Rxn_module_mGLuR.xml",
                              "Rxn_module_IP3R.xml"]



if __name__ == "__main__":
    
    # 1 no mGluR no RyR
   

    my_rxn_f = read_in_files(flist_ER_no_IP3)
    with  open("Rxn_ER_no_IP3.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))

    my_rxn_f = read_in_files(flist_ER)
    with  open("Rxn_start.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
