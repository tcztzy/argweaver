use nom::{
    bytes::complete::tag,
    character::complete::{alphanumeric1, newline, tab},
    multi::separated_list1,
    sequence::{delimited, terminated},
    IResult,
};

pub fn parse_names<'a>(input: &'a str) -> IResult<&'a str, Vec<&'a str>> {
    delimited(
        terminated(tag("NAMES"), tab),
        separated_list1(tab, alphanumeric1),
        newline,
    )(input)
}

#[test]
fn test_parse_names() {
    let input = "NAMES\tA\tB\tC\n";
    let res = parse_names(input);
    assert_eq!(res, Ok(("", vec!["A", "B", "C"])));
}
