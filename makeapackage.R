#first, establish the connection to github with keys
#on xenon4:
git config --global user.name xxxmichixxx
git config --global user.email michaela.c.schwaiger@gmail.com
git config --list
ssh-keygen
less .ssh/id_rsa.pub #copy this public key to github (in my profile settings)
ssh -T git@github.com #this tests the connection

#generate github repository on github, and copy ssh url
#start a new project in rstudio and copy in github key
#add .gitignore and MiniChip.Rproj to .gitignore file (so it doesnt go to github)

#to create an R package skeleton:
usethis::create_package("MiniChip")
usethis::use_roxygen_md()
#in rstudio, under Project Options, Buid Tools, check boxes to make sure Roxygen is used
usethis::use_news_md()
usethis::use_testthat()


# make a file functionname.R for each function in the R directory, add @blabla stuff for roxygen
#make a file test-functionname.R and add tests for this function
#edit DESCRIPTION file:
# import packages, do not add them under depends
usethis::use_package("ComplexHeatmap") #this will automatically add it to imports in description file

#when edits are made, go to Build, click: Install and Restart
#then add the documentations:
devtools::document() # this crashesR
#better use the document tab (just below test tab)
#then install and restart again
#to test package:more,TestPackage
#to really test package (slower): Check
# to deliver package to users: Build, Build Source package: this makes a tar file which can be sent by email and installed with install.packages() by the user

#to add to git repository: Commit
#to add to Github: Push

#add a vignette
usethis::use_vignette("MiniChip-vignette")

