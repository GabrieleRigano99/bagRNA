Submit-block ::= {
  contact {
    contact {
      name name {
        last "Rigano",
        first "Gabriele"
      },
      affil std {
        affil "University of Rome Tor Vergata",
        div "Department of Biology",
        city "Rome",
        sub "University of Rome Tor Vergata",
        country "Italy",
        street "Via della ricerca scientifica",
        email "gabriele.rigano@unitn.it",
        postal-code "00133"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Rigano",
            first "Gabriele",
            initials "G.R.",
            suffix ""
          }
        },
        {
          name name {
            last "Bonomo",
            first "Andrea",
            initials "A.B.",
            suffix ""
          }
        }
      },
      affil std {
        affil "University of Rome Tor Vergata",
        div "Department of Biology",
        city "Rome",
        sub "University of Rome Tor Vergata",
        country "Italy",
        street "Via della ricerca scientifica",
        postal-code "00133"
      }
    }
  },
  subtype new
}

Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Rigano",
              first "Gabriele",
              initials "G.R.",
              suffix ""
            }
          },
          {
            name name {
              last "Bonomo",
              first "Andrea",
              initials "A.B.",
              suffix ""
            }
          }
        }
      },
      title "Genome annotation of Sporothrix schenckii 1099-18"
    }
  }
}

Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        ""
      }
    }
  }
}

