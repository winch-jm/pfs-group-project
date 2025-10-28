package main
import(
	"fmt"
)
func main(){
	cells := ReadCSV("data/ctl_subset.csv")
	fmt.Print(cells[0].name)

}