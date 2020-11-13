import argparse

from idiva.clf.utils import create_df, get_clf

from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf_file', type=Path, help="VCF file to classify the variants")
    parser.add_argument('--classifier', type=str, default='', help="Classifier to use for predictions")
    flags = parser.parse_args()

    assert isinstance(flags.vcf_file, Path)
    assert flags.vcf_file.is_file()

    flags.classifier
    # TODO: check validity

    return dict(vcf=flags.vcf_file, classifier=flags.classifier)


def predict(vcf: Path, classifier: str = None):
    """
    Reads a given vcf file, transforms it into a dataframe and uses a pretrained classifier to classify each variant in
    the vcf with 0 for control and 1 for "sick". It then adds the label to the input vcf and returns it (or saves it).
    """

    assert vcf.is_file()

    df = create_df(vcf)
    clf = get_clf(classifier)
    # TODO: Don't bother with numpy arrays unless there is a good reason -- ?
    prediction = clf.predict(df)

    return prediction


def main():
    prediction = predict(**parse_args())
    # todo: somehow include the predictions into the input vcf file


if __name__ == '__main__':
    main()
