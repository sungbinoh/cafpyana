POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      EXTENSION="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      OUTPUT="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev.sh
. /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup root v6_28_12 -q e26:p3915:prof
root 'MakesBruce.C("test.root", "output.root")
