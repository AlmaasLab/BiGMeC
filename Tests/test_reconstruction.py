# test_reconstruction
import sys
sys.path.append("..\\Scripts")
import bigmec

def test_bigmec():
    summary_df = bigmec.run("../Data/mibig/1.gbk", '../Data/constructed_pathways/', '../Models/BiGG_universal_model.json')
    assert summary_df["Success"][0] == 1

if __name__ == '__main__':
    test_bigmec()
